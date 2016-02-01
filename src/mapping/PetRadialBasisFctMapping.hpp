#pragma once
#ifndef PRECICE_NO_PETSC

#include "mapping/Mapping.hpp"
#include "impl/BasisFunctions.hpp"
#include "tarch/la/DynamicVector.h"
#include "utils/MasterSlave.hpp"
#include "utils/Petsc.hpp"
#include <limits>
#include <typeinfo>

#include "petnum.hpp"
#include "petscmat.h"
#include "petscksp.h"
#include "petsclog.h"

#include <iostream>
using std::cout;
using std::endl;
#include "utils/prettyprint.hpp"
#include "utils/EventTimings.hpp"


namespace precice {
namespace mapping {

/**
 * @brief Mapping with radial basis functions using the Petsc library to solve the resulting system.
 *
 * With help of the input data points and values an interpolant is constructed.
 * The interpolant is formed by a weighted sum of conditionally positive radial
 * basis functions and a (low order) polynomial, and evaluated at the output
 * data points.
 *
 * The radial basis function type has to be given as template parameter, and has
 * to be one of the defined types in this file.
 */
template<typename RADIAL_BASIS_FUNCTION_T>
class PetRadialBasisFctMapping : public Mapping
{
public:

  /**
   * @brief Constructor.
   *
   * @param[in] constraint Specifies mapping to be consistent or conservative.
   * @param[in] function Radial basis function used for mapping.
   * @param[in] solverRto Relative tolerance for the linear solver.
   *
   * For description on convergence testing and meaning of solverRtol see http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPConvergedDefault.html#KSPConvergedDefault
   */
  PetRadialBasisFctMapping (
    Constraint              constraint,
    int                     dimensions,
    RADIAL_BASIS_FUNCTION_T function,
    bool                    xDead,
    bool                    yDead,
    bool                    zDead,
    double                  solverRtol = 1e-9);

  /// Destroys the Petsc KSP and the _deadAxis array
  virtual ~PetRadialBasisFctMapping();

  /// @brief Computes the mapping coefficients from the in- and output mesh.
  virtual void computeMapping();

  /// @brief Returns true, if computeMapping() has been called.
  virtual bool hasComputedMapping() const;

  /// @brief Removes a computed mapping.
  virtual void clear();

  /// @brief Maps input data to output data from input mesh to output mesh.
  virtual void map(int inputDataID, int outputDataID);

private:

  /// @brief Logging device.
  static tarch::logging::Log _log;

  bool _hasComputedMapping;

  /// @brief Radial basis function type used in interpolation.
  RADIAL_BASIS_FUNCTION_T _basisFunction;

  petsc::Matrix _matrixC;

  petsc::Matrix _matrixA;

  KSP _solver;

  ISLocalToGlobalMapping _ISmapping;

  double _solverRtol;

  /// true if the mapping along some axis should be ignored
  bool* _deadAxis;

  /// Deletes all dead directions from fullVector and returns a vector of reduced dimensionality.
  // utils::DynVector reduceVector(const utils::DynVector& fullVector);

  virtual bool doesVertexContribute(int vertexID) const;
};

// --------------------------------------------------- HEADER IMPLEMENTATIONS

template<typename RADIAL_BASIS_FUNCTION_T>
tarch::logging::Log PetRadialBasisFctMapping<RADIAL_BASIS_FUNCTION_T>::_log("precice::mapping::PetRadialBasisFctMapping");

template<typename RADIAL_BASIS_FUNCTION_T>
PetRadialBasisFctMapping<RADIAL_BASIS_FUNCTION_T>::PetRadialBasisFctMapping
(
  Constraint              constraint,
  int                     dimensions,
  RADIAL_BASIS_FUNCTION_T function,
  bool                    xDead,
  bool                    yDead,
  bool                    zDead,
  double                  solverRtol)
  :
  Mapping ( constraint, dimensions ),
  _hasComputedMapping ( false ),
  _basisFunction ( function ),
  _matrixC(PETSC_COMM_WORLD, "C"),
  _matrixA(PETSC_COMM_WORLD, "A"),
  _solverRtol(solverRtol)
{
  setInputRequirement(VERTEX);
  setOutputRequirement(VERTEX);
  _deadAxis = new bool[dimensions];

  if (getDimensions()==2) {
    _deadAxis[0] = xDead;
    _deadAxis[1] = yDead;
    preciceCheck(not (xDead && yDead),
                 "setDeadAxis()", "You cannot choose all axis to be dead for a RBF mapping");
    preciceCheck(not zDead,
                 "setDeadAxis()", "You cannot dead out the z axis if dimension is set to 2");
  }
  else if (getDimensions()==3) {
    _deadAxis[0] = xDead;
    _deadAxis[1] = yDead;
    _deadAxis[2] = zDead;
    preciceCheck(not (xDead && yDead && zDead), "setDeadAxis()", "You cannot  "
                 << " choose all axis to be dead for a RBF mapping");
  }
  else {
    assertion(false);
  }

  KSPCreate(PETSC_COMM_WORLD, &_solver);
}

template<typename RADIAL_BASIS_FUNCTION_T>
PetRadialBasisFctMapping<RADIAL_BASIS_FUNCTION_T>::~PetRadialBasisFctMapping()
{
  delete[] _deadAxis;
  // PetscErrorCode ierr = 0;
  // Commenting out the next line most likely introduces a memory leak
  // However, not commenting it out introduces a memory error, which remains untraceable
  // ierr = KSPDestroy(&_solver); CHKERRV(ierr);
}



template<typename RADIAL_BASIS_FUNCTION_T>
void PetRadialBasisFctMapping<RADIAL_BASIS_FUNCTION_T>::computeMapping()
{
  preciceTrace("computeMapping()");
  precice::utils::Event e(__func__);

  assertion2(input()->getDimensions() == output()->getDimensions(),
             input()->getDimensions(), output()->getDimensions());
  int dimensions = input()->getDimensions();
  mesh::PtrMesh inMesh;
  mesh::PtrMesh outMesh;
  if (getConstraint() == CONSERVATIVE){
    inMesh = output();
    outMesh = input();
  }
  else {
    inMesh = input();
    outMesh = output();
  }

  int deadDimensions = 0;
  for (int d=0; d<dimensions; d++) {
    if (_deadAxis[d]) deadDimensions +=1;
  }
  size_t polyparams = 1 + dimensions - deadDimensions;

  // Indizes that are used to build the Petsc Index set
  std::vector<int> myIndizes;

  // Indizes for Q^T, holding the polynom
  if (utils::MasterSlave::_rank <= 0) // Rank 0 or not in MasterSlave mode
    for (size_t i = 0; i < polyparams; i++)
      myIndizes.push_back(i); // polyparams reside in the first rows (which are always on rank 0)

  // Indizes for the vertices with polyparams offset
  for (const mesh::Vertex& v : inMesh->vertices())
    if (v.isOwner())
      myIndizes.push_back(v.getGlobalIndex() + polyparams);

  auto inputSize = myIndizes.size();
  auto n = inputSize; // polyparams, if on rank 0, are included here

  if (utils::MasterSlave::_rank <= 0)
    inputSize -= polyparams; // Subtract polyparams on rank 0, so we only have number of vertices.

  auto outputSize = outMesh->vertices().size();

  PetscErrorCode ierr = 0;
  preciceDebug("inMesh->vertices().size() = " << inMesh->vertices().size());
  preciceDebug("outMesh->vertices().size() = " << outMesh->vertices().size());

  // Create a symmetric, sparse matrix with n x n local size.
  _matrixC.reset();
  _matrixC.init(n, n, PETSC_DETERMINE, PETSC_DETERMINE, MATSBAIJ);
  preciceDebug("Set matrix C to local size " << n << " x " << n);
  ierr = MatSetOption(_matrixC.matrix, MAT_SYMMETRY_ETERNAL, PETSC_TRUE); CHKERRV(ierr);

  // Create a matrix with outputSize x n local size.
  _matrixA.reset();
  _matrixA.init(outputSize, n, PETSC_DETERMINE, PETSC_DETERMINE, MATAIJ);
  preciceDebug("Set matrix A to local size " << outputSize << " x " << n);

  KSPReset(_solver);

  IS ISlocal, ISglobal, ISidentity, ISidentityGlobal;
  ISLocalToGlobalMapping ISidentityMapping;
  // Create an index set which maps myIndizes to continous chunks of matrix rows.
  ierr = ISCreateGeneral(PETSC_COMM_WORLD, myIndizes.size(), myIndizes.data(), PETSC_COPY_VALUES, &ISlocal); CHKERRV(ierr);
  ierr = ISAllGather(ISlocal, &ISglobal); CHKERRV(ierr); // Gather the IS from all processors
  ierr = ISLocalToGlobalMappingCreateIS(ISglobal, &_ISmapping); CHKERRV(ierr); // Make it a mapping
  
  // Create an identity mapping and use that for the rows of matrixA.
  ierr = ISCreateStride(PETSC_COMM_WORLD, _matrixA.ownerRange().second - _matrixA.ownerRange().first, _matrixA.ownerRange().first, 1, &ISidentity); CHKERRV(ierr);
  ISSetIdentity(ISidentity);
  ierr = ISAllGather(ISidentity, &ISidentityGlobal); CHKERRV(ierr);
  ierr = ISLocalToGlobalMappingCreateIS(ISidentityGlobal, &ISidentityMapping); CHKERRV(ierr);

  ierr = MatSetLocalToGlobalMapping(_matrixC.matrix, _ISmapping, _ISmapping); CHKERRV(ierr); // Set mapping for rows and cols
  ierr = MatSetLocalToGlobalMapping(_matrixA.matrix, ISidentityMapping, _ISmapping); CHKERRV(ierr); // Set mapping only for cols, use identity for rows

  // if (utils::MasterSlave::_rank <= 0) {
  //   cout << "== ISidentityGlobal ==" << endl;
  //   ISView(ISidentityGlobal, PETSC_VIEWER_STDOUT_SELF);
  //   cout << "== ISglobal ==" << endl;
  //   ISView(ISglobal, PETSC_VIEWER_STDOUT_SELF);
  // }

  int i = 0;
  utils::DynVector distance(dimensions);

  // We do preallocating of the matrices C and A. That means we traverse the input data once, just
  // to know where we have entries in the sparse matrix. This information petsc can use to
  // preallocate the matrix. In the second phase we actually fill the matrix.

  // -- BEGIN PREALLOC LOOP FOR MATRIX C --
  // TODO Testen ob Preallocation perfekt ist.
  // int logPreallocCLoop = 1;
  // PetscLogEventRegister("Prealloc Matrix C", 0, &logPreallocCLoop);
  // PetscLogEventBegin(logPreallocCLoop, 0, 0, 0, 0);
  // PetscInt nnz[n]; // Number of non-zeros per row
  // unsigned int totalNNZ = 0; // Total number of non-zeros in matrix
  // for (const mesh::Vertex& iVertex : inMesh->vertices()) {
  //   int currentRow = iVertex.getGlobalIndex();
  //   nnz[currentRow] = 0;
  //   for (int j=iVertex.getID(); j < inputSize; j++) {
  //     if (currentRow == iVertex.getGlobalIndex()) {
  //       nnz[currentRow]++; // Since we need to set at least zeros on the main diagonal, we have an entry there. No further test necessary.
  //       continue;
  //     }
  //     distance = iVertex.getCoords() - inMesh->vertices()[j].getCoords();
  //     double coeff = _basisFunction.evaluate(norm2(distance));
  //     if ( not equals(coeff, 0.0)) {
  //       nnz[currentRow]++;
  //     }
  //   }
  //   nnz[currentRow] += dimensions + 1 - deadDimensions; // coefficients of the linear polynom
  //   totalNNZ += nnz[currentRow];
  //   i++;
  // }
  // for (int r = inputSize; r < n; r++) {
  //   nnz[r] = 1; // iterate through the lower part, no coefficients there, but 0 on the main diagonal
  // }
  // PetscLogEventEnd(logPreallocCLoop, 0, 0, 0, 0);
  // i = 0;
  // // -- END PREALLOC LOOP FOR MATRIX C --

  // ierr = MatSeqSBAIJSetPreallocation(_matrixC.matrix, 1, PETSC_DEFAULT, nnz);

  // -- BEGIN FILL LOOP FOR MATRIX C --
  int logCLoop = 2;
  PetscLogEventRegister("Filling Matrix C", 0, &logCLoop);
  PetscLogEventBegin(logCLoop, 0, 0, 0, 0);
  // We collect entries for each row and set them blockwise using MatSetValues.
  for (const mesh::Vertex& inVertex : inMesh->vertices()) {
    if (not inVertex.isOwner())
      continue;

    int row = inVertex.getGlobalIndex() + polyparams;

    // -- SETS THE POLYNOM PART OF THE MATRIX --
    PetscInt colIdx[_matrixC.getSize().second];     // holds the columns indices of the entries
    PetscScalar colVals[_matrixC.getSize().second]; // holds the values of the entries
    PetscInt polyRow = 0, polyCol = row;
    PetscScalar y = 1;
    ierr = MatSetValuesLocal(_matrixC.matrix, 1, &polyRow, 1, &polyCol, &y, INSERT_VALUES); CHKERRV(ierr);
    int actualDim = 0;
    for (int dim = 0; dim < dimensions; dim++) {
      if (not _deadAxis[dim]) {
        y = inVertex.getCoords()[dim];
        polyRow = 1 + actualDim;
        actualDim++;
        ierr = MatSetValuesLocal(_matrixC.matrix, 1, &polyRow, 1, &polyCol, &y, INSERT_VALUES); CHKERRV(ierr);
      }
    }

    // -- SETS THE COEFFICIENTS --
    PetscInt colNum = 0;  // holds the number of columns
    for (mesh::Vertex& vj : inMesh->vertices()) {
      distance = inVertex.getCoords() - vj.getCoords();
      for (int d = 0; d < dimensions; d++) {
        if (_deadAxis[d]) {
          distance[d] = 0;
        }
      }
      double coeff = _basisFunction.evaluate(norm2(distance));
      if (not tarch::la::equals(coeff, 0.0)) {
        colVals[colNum] = coeff;
        colIdx[colNum] = vj.getGlobalIndex() + polyparams; // column of entry is the globalIndex
        colNum++;
      }
      #ifdef Asserts
      if (coeff == std::numeric_limits<double>::infinity()) {
        preciceError("computeMapping()", "C matrix element has value inf. "
                     << "i = " << i
                     << ", coords i = " << inVertex.getCoords() << ", coords j = "
                     << vj.getCoords() << ", dist = "
                     << distance << ", norm2 = " << norm2(distance) << ", rbf = "
                     << coeff
                     << ", rbf type = " << typeid(_basisFunction).name());
      }
      # endif
    }

    ierr = MatSetValuesLocal(_matrixC.matrix, 1, &row, colNum, colIdx, colVals, INSERT_VALUES); CHKERRV(ierr);
    i++;
  }
  PetscLogEventEnd(logCLoop, 0, 0, 0, 0);
  // -- END FILL LOOP FOR MATRIX C --

  // Petsc requires that all diagonal entries are set, even if set to zero.
  _matrixC.assemble(MAT_FLUSH_ASSEMBLY);
  petsc::Vector zeros(_matrixC);
  MatDiagonalSet(_matrixC.matrix, zeros.vector, ADD_VALUES);

  // Begin assembly here, all assembly is ended at the end of this function.
  ierr = MatAssemblyBegin(_matrixC.matrix, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);

  // -- BEGIN PREALLOC LOOP FOR MATRIX A --
  // int logPreallocALoop = 3;
  // PetscLogEventRegister("Prealloc Matrix A", 0, &logPreallocALoop);
  // PetscLogEventBegin(logPreallocALoop, 0, 0, 0, 0);
  // PetscInt nnzA[outputSize];
  // i = 0;
  // for (const mesh::Vertex& iVertex : outMesh->vertices()) {
  //   nnzA[i] = 0;
  //   for (const mesh::Vertex& jVertex : inMesh->vertices()) {
  //     distance = iVertex.getCoords() - jVertex.getCoords();
  //     double coeff = _basisFunction.evaluate(norm2(distance));
  //     if ( not equals(coeff, 0.0)) {
  //       nnzA[i]++;
  //     }
  //   }
  //   nnzA[i] += dimensions + 1 - deadDimensions;
  //   i++;
  // }
  // PetscLogEventEnd(logPreallocALoop, 0, 0, 0, 0);
  // // -- END PREALLOC LOOP FOR MATRIX A --

  // ierr = MatSeqAIJSetPreallocation(_matrixA.matrix, PETSC_DEFAULT, nnzA); CHKERRV(ierr);

  // -- BEGIN FILL LOOP FOR MATRIX A --
  int logALoop = 4;
  PetscLogEventRegister("Filling Matrix A", 0, &logALoop);
  PetscLogEventBegin(logALoop, 0, 0, 0, 0);
  i = 0;
  for (int it = _matrixA.ownerRange().first; it < _matrixA.ownerRange().second; it++) {
    // hier colIdx, colVals über inMesh->vertices.count dimensionieren?
    PetscInt colNum = 0;
    PetscInt colIdx[_matrixC.getSize().second];     // holds the columns indices of the entries
    PetscScalar colVals[_matrixC.getSize().second]; // holds the values of the entries
    const mesh::Vertex& oVertex = outMesh->vertices()[it - _matrixA.ownerRange().first];
    // preciceDebug("Matrix A, Row = " << it);

    // -- SET THE POLYNOM PART OF THE MATRIX --
    PetscInt polyRow = it, polyCol = 0;
    PetscScalar y = 1;
    ierr = MatSetValuesLocal(_matrixA.matrix, 1, &polyRow, 1, &polyCol, &y, INSERT_VALUES); CHKERRV(ierr);
    for (int dim = 0; dim < dimensions; dim++) {
      if (not _deadAxis[dim]) {
        y = oVertex.getCoords()[dim];
        polyCol++;
        ierr = MatSetValuesLocal(_matrixA.matrix, 1, &polyRow, 1, &polyCol, &y, INSERT_VALUES); CHKERRV(ierr);
      }
    }

    // -- SETS THE COEFFICIENTS --
    int j = 0;
    for (const mesh::Vertex& inVertex : inMesh->vertices()) {
      distance = oVertex.getCoords() - inVertex.getCoords();
      for (int d = 0; d < dimensions; d++) {
        if (_deadAxis[d])
          distance[d] = 0;
      }
      double coeff = _basisFunction.evaluate(norm2(distance));
      if ( not tarch::la::equals(coeff, 0.0)) {
        colVals[colNum] = coeff;
        colIdx[colNum] = inVertex.getGlobalIndex() + polyparams;
        colNum++;
      }

      #     ifdef Asserts
      if (coeff == std::numeric_limits<double>::infinity()){
        preciceError("computeMapping()", "A matrix element has value inf. "
                     << "i = " << i << ", j = " << j
                     << ", coords i = " << oVertex.getCoords() << ", coords j = "
                     << inVertex.getCoords() << ", dist = "
                     << distance << ", norm2 = " << norm2(distance) << ", rbf = "
                     << coeff
                     << ", rbf type = " << typeid(_basisFunction).name());
      }
      #     endif
      j++;
    }
    ierr = MatSetValuesLocal(_matrixA.matrix, 1, &it, colNum, colIdx, colVals, INSERT_VALUES); CHKERRV(ierr);
    i++;
  }
  PetscLogEventEnd(logALoop, 0, 0, 0, 0);
  // -- END FILL LOOP FOR MATRIX A --

  ierr = MatAssemblyBegin(_matrixA.matrix, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
  ierr = MatAssemblyEnd(_matrixC.matrix, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
  ierr = MatAssemblyEnd(_matrixA.matrix, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
  KSPSetOperators(_solver, _matrixC.matrix, _matrixC.matrix);
  KSPSetTolerances(_solver, _solverRtol, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
  KSPSetFromOptions(_solver);

  // if (totalNNZ > static_cast<size_t>(20*n)) {
  //   preciceDebug("Using Cholesky decomposition as direct solver for dense matrix.");
  //   PC prec;
  //   KSPSetType(_solver, KSPPREONLY);
  //   KSPGetPC(_solver, &prec);
  //   PCSetType(prec, PCCHOLESKY);
  //   PCFactorSetShiftType(prec, MAT_SHIFT_NONZERO);
  // }

  _hasComputedMapping = true;
  // _matrixA.view();
}

template<typename RADIAL_BASIS_FUNCTION_T>
bool PetRadialBasisFctMapping<RADIAL_BASIS_FUNCTION_T>:: hasComputedMapping() const
{
  return _hasComputedMapping;
}

template<typename RADIAL_BASIS_FUNCTION_T>
void PetRadialBasisFctMapping<RADIAL_BASIS_FUNCTION_T>:: clear()
{
  preciceTrace("clear()");
  _matrixC.reset();
  _matrixA.reset();
  _hasComputedMapping = false;
  // ISmapping destroy??
}

template<typename RADIAL_BASIS_FUNCTION_T>
void PetRadialBasisFctMapping<RADIAL_BASIS_FUNCTION_T>:: map
(
  int inputDataID,
  int outputDataID )
{
  precice::utils::Event e(__func__);
  preciceTrace2("map()", inputDataID, outputDataID);

  assertion(_hasComputedMapping);
  assertion2(input()->getDimensions() == output()->getDimensions(),
             input()->getDimensions(), output()->getDimensions());
  using namespace tarch::la;
  PetscErrorCode ierr = 0;
  KSPConvergedReason convReason;
  utils::DynVector& inValues = input()->data(inputDataID)->values();
  utils::DynVector& outValues = output()->data(outputDataID)->values();

  int valueDim = input()->data(inputDataID)->getDimensions();
  assertion2(valueDim == output()->data(outputDataID)->getDimensions(),
             valueDim, output()->data(outputDataID)->getDimensions());
  int deadDimensions = 0;
  for (int d=0; d<getDimensions(); d++) {
    if (_deadAxis[d]) deadDimensions +=1;
  }
  int polyparams = 1 + getDimensions() - deadDimensions;

  if (getConstraint() == CONSERVATIVE) {
    preciceDebug("Map conservative");
    static int mappingIndex = 0;
    petsc::Vector Au(_matrixC, "Au");
    petsc::Vector out(_matrixC, "out");
    petsc::Vector in(_matrixA, "in");
    ierr = VecSetLocalToGlobalMapping(in.vector, _ISmapping); CHKERRV(ierr);

    for (int dim=0; dim < valueDim; dim++) {
      preciceDebug("in vector ownerRange = " << in.ownerRange());
      for (int i = in.ownerRange().first; i < in.ownerRange().second; i++) {
        int index = i - in.ownerRange().first; // Relative (local) index
        int globalIndex = input()->vertices()[index].getGlobalIndex();
        preciceDebug("Filling input vector(" << globalIndex << ") = " <<  inValues[index*valueDim + dim]);
        VecSetValueLocal(in.vector, globalIndex, inValues[index*valueDim + dim], INSERT_VALUES);        // Dies besser als VecSetValuesLocal machen
      }
      in.assemble();
      // in.view();
      ierr = MatMultTranspose(_matrixA.matrix, in.vector, Au.vector); CHKERRV(ierr);
      ierr = KSPSolve(_solver, Au.vector, out.vector); CHKERRV(ierr);
      ierr = KSPGetConvergedReason(_solver, &convReason); CHKERRV(ierr);
      if (convReason < 0) {
        preciceError(__func__, "RBF linear system has not converged.");
      }
      VecChop(out.vector, 1e-9);
      // Copy mapped data to output data values
      const PetscScalar *outArray;
      ierr = VecGetArrayRead(out.vector, &outArray);
      int size = out.getLocalSize();
      preciceDebug("Local out vector size = " << size);
      for (int i=0; i < size-polyparams; i++){
        outValues[i*valueDim + dim] = outArray[i+polyparams]; // hier noch das index set beachten?
      }
      VecRestoreArrayRead(out.vector, &outArray);
    }
    mappingIndex++;
  }
  else { // Map consistent
    preciceDebug("Map consistent");
    petsc::Vector p(_matrixC, "p");
    petsc::Vector in(_matrixC, "in");
    petsc::Vector out(_matrixA, "out");
    ierr = VecSetLocalToGlobalMapping(in.vector, _ISmapping); CHKERRV(ierr);
    const PetscScalar *vecArray;

    // For every data dimension, perform mapping
    for (int dim=0; dim < valueDim; dim++){
      // Fill input from input data values
      preciceDebug("in vector ownerRange = " << in.ownerRange());
      preciceDebug("polyparams = " << polyparams << ", valueDim = " << valueDim << ", dim = " << dim);
      // for (int i = in.ownerRange().first; i < in.ownerRange().second; i++) {
      auto localPolyparams = utils::MasterSlave::_rank > 0 ? 0 : polyparams; // Set localPolyparams only when root rank
      for (int i = in.ownerRange().first + localPolyparams; i < in.ownerRange().second; i++) {
        preciceDebug("Begin Loop, i = " << i);

        // if (i < polyparams) // The polyparams remain zero, skipping.
          // continue;
        // int index = i - in.ownerRange().first; // Relative (local) index
        // preciceDebug("local index = " << index);
        int globalIndex = input()->vertices()[i-polyparams].getGlobalIndex(); // i - ownerRange.first ?
        preciceDebug("globalIndex = " << globalIndex);
        preciceDebug("Filling input vector(" << globalIndex+polyparams << ") = inValues[" << (i-polyparams)*valueDim + dim << "] = " << inValues[(i-polyparams)*valueDim + dim]);
        VecSetValueLocal(in.vector, globalIndex+polyparams, inValues[(i-polyparams)*valueDim + dim], INSERT_VALUES);        // Dies besser als VecSetValuesLocal machen

        // Begin Benjamin
        // preciceDebug("Filling input vector(" << i << ") = " <<  inValues[(i-polyparams)*valueDim + dim]);
        // VecSetValueLocal(in.vector, i, inValues[(i-polyparams)*valueDim + dim], INSERT_VALUES);        // Dies besser als VecSetValuesLocal machen
        // End Benjamin
        preciceDebug("End Loop, i = " << i);

      }
      preciceDebug("Finished in vector construction.")
      in.assemble();
      // in.view();
      ierr = KSPSolve(_solver, in.vector, p.vector); CHKERRV(ierr);
      ierr = KSPGetConvergedReason(_solver, &convReason); CHKERRV(ierr);
      if (convReason < 0) {
        preciceError(__func__, "RBF linear system has not converged.");
      }
      ierr = MatMult(_matrixA.matrix, p.vector, out.vector); CHKERRV(ierr);
      VecChop(out.vector, 1e-9);
      // Copy mapped data to output data values
      ierr = VecGetArrayRead(out.vector, &vecArray);
      int size = out.getLocalSize();
      preciceDebug("Local out vector size = " << size);
      for (int i=0; i < size; i++) {
        outValues[i*valueDim + dim] = vecArray[i];
      }
      VecRestoreArrayRead(out.vector, &vecArray);
    }
  }
}


template<typename RADIAL_BASIS_FUNCTION_T>
bool PetRadialBasisFctMapping<RADIAL_BASIS_FUNCTION_T>::doesVertexContribute(int vertexID) const
{
  // FIXME: Use a sane calculation here
  preciceTrace(__func__);
  preciceDebug("Mesh Size = " << output()->vertices().size());
  // for (auto& v : output()->vertices())
  // if (not v.isOwner())
  // preciceDebug("!!! NOT OWNER !!!");
  // if (v.getID() == vertexID)
  // preciceDebug("Owner " << v.isOwner());
  if (not _basisFunction.hasCompactSupport())
    return true;

  return true;

  
}


}} // namespace precice, mapping

#endif
