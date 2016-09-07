#pragma once
#ifndef PRECICE_NO_PETSC

#include <limits>
#include <typeinfo>
#include <map>

#include "mapping/Mapping.hpp"
#include "impl/BasisFunctions.hpp"
#include "tests/PetRadialBasisFctMappingTest.hpp"
#include "utils/MasterSlave.hpp"
#include "utils/Petsc.hpp"
namespace petsc = precice::utils::petsc;

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
 * The radial basis function type has to be given as template parameter.
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
   * @param[in] solverRtol Relative tolerance for the linear solver.
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

  /// Deletes the PETSc objects and the _deadAxis array
  virtual ~PetRadialBasisFctMapping();

  /// Computes the mapping coefficients from the in- and output mesh.
  virtual void computeMapping() override;

  /// Returns true, if computeMapping() has been called.
  virtual bool hasComputedMapping() const override;

  /// Removes a computed mapping.
  virtual void clear() override;

  /// Maps input data to output data from input mesh to output mesh.
  virtual void map(int inputDataID, int outputDataID) override;

  friend class precice::mapping::tests::PetRadialBasisFctMappingTest;

private:

  /// Logging device.
  static logging::Logger _log;

  bool _hasComputedMapping;

  /// Radial basis function type used in interpolation.
  RADIAL_BASIS_FUNCTION_T _basisFunction;

  petsc::Matrix _matrixC;

  petsc::Matrix _matrixA;

  KSP _solver;

  ISLocalToGlobalMapping _ISmapping;

  double _solverRtol;

  /// true if the mapping along some axis should be ignored
  bool* _deadAxis;

  virtual bool doesVertexContribute(int vertexID) const override;

  void incPrealloc(PetscInt* diag, PetscInt* offDiag, int pos, int begin, int end);

  /// Stores the solution from the previous iteration
  std::map<unsigned int, petsc::Vector> previousSolution;
};

// --------------------------------------------------- HEADER IMPLEMENTATIONS

template<typename RADIAL_BASIS_FUNCTION_T>
logging::Logger PetRadialBasisFctMapping<RADIAL_BASIS_FUNCTION_T>::_log("precice::mapping::PetRadialBasisFctMapping");

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
  PetscBool petscIsInitialized;
  PetscErrorCode ierr = 0;
  PetscInitialized(&petscIsInitialized);
  if (petscIsInitialized) {
    ierr = ISLocalToGlobalMappingDestroy(&_ISmapping); CHKERRV(ierr);
    ierr = KSPDestroy(&_solver); CHKERRV(ierr);
  }
}



template<typename RADIAL_BASIS_FUNCTION_T>
void PetRadialBasisFctMapping<RADIAL_BASIS_FUNCTION_T>::computeMapping()
{
  preciceTrace("computeMapping()");
  precice::utils::Event e(__func__);

  assertion(input()->getDimensions() == output()->getDimensions(),
             input()->getDimensions(), output()->getDimensions());
  int dimensions = input()->getDimensions();
  mesh::PtrMesh inMesh;
  mesh::PtrMesh outMesh;
  if (getConstraint() == CONSERVATIVE) {
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

  // Matrix C: Symmetric, sparse matrix with n x n local size.
  _matrixC.reset();
  _matrixC.init(n, n, PETSC_DETERMINE, PETSC_DETERMINE, MATSBAIJ);
  preciceDebug("Set matrix C to local size " << n << " x " << n);
  ierr = MatSetOption(_matrixC, MAT_SYMMETRY_ETERNAL, PETSC_TRUE); CHKERRV(ierr);

  // Matrix A: Sparse matrix with outputSize x n local size.
  _matrixA.reset();
  _matrixA.init(outputSize, n, PETSC_DETERMINE, PETSC_DETERMINE, MATAIJ);
  preciceDebug("Set matrix A to local size " << outputSize << " x " << n);

  KSPReset(_solver);

  const int ownerRangeABegin = _matrixA.ownerRange().first;
  const int ownerRangeAEnd = _matrixA.ownerRange().second;
  
  IS ISlocal, ISlocalInv, ISglobal, ISidentity, ISidentityGlobal;
  ISLocalToGlobalMapping ISidentityMapping;
  // Create an index set which maps myIndizes to continous chunks of matrix rows.
  ierr = ISCreateGeneral(PETSC_COMM_WORLD, myIndizes.size(), myIndizes.data(), PETSC_COPY_VALUES, &ISlocal); CHKERRV(ierr);
  ierr = ISSetPermutation(ISlocal); CHKERRV(ierr);
  ierr = ISInvertPermutation(ISlocal, myIndizes.size(), &ISlocalInv); CHKERRV(ierr);
  ierr = ISAllGather(ISlocalInv, &ISglobal); CHKERRV(ierr); // Gather the IS from all processors
  ierr = ISLocalToGlobalMappingCreateIS(ISglobal, &_ISmapping); CHKERRV(ierr); // Make it a mapping
  
  // Create an identity mapping and use that for the rows of matrixA.
  ierr = ISCreateStride(PETSC_COMM_WORLD, ownerRangeAEnd - ownerRangeABegin, ownerRangeABegin, 1, &ISidentity); CHKERRV(ierr);
  ierr = ISSetIdentity(ISidentity); CHKERRV(ierr);
  ierr = ISAllGather(ISidentity, &ISidentityGlobal); CHKERRV(ierr);
  ierr = ISLocalToGlobalMappingCreateIS(ISidentityGlobal, &ISidentityMapping); CHKERRV(ierr);

  ierr = MatSetLocalToGlobalMapping(_matrixC, _ISmapping, _ISmapping); CHKERRV(ierr); // Set mapping for rows and cols
  ierr = MatSetLocalToGlobalMapping(_matrixA, ISidentityMapping, _ISmapping); CHKERRV(ierr); // Set mapping only for cols, use identity for rows

  // Destroy all local index sets and mappings
  ierr = ISDestroy(&ISlocal); CHKERRV(ierr);
  ierr = ISDestroy(&ISlocalInv); CHKERRV(ierr);
  ierr = ISDestroy(&ISglobal); CHKERRV(ierr);
  ierr = ISDestroy(&ISidentity); CHKERRV(ierr);
  ierr = ISDestroy(&ISidentityGlobal); CHKERRV(ierr);
  ierr = ISLocalToGlobalMappingDestroy(&ISidentityMapping); CHKERRV(ierr);

  utils::DynVector distance(dimensions);

  // We do preallocating of the matrices C and A. That means we traverse the input data once, just
  // to know where we have entries in the sparse matrix. This information petsc can use to
  // preallocate the matrix. In the second phase we actually fill the matrix.

  // -- BEGIN PREALLOC LOOP FOR MATRIX C --
  // TODO Testen ob Preallocation perfekt ist.
  // MatSetOption(_matrixC, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  // MatSetOption(_matrixA, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  // {
  //   preciceDebug("Begin preallocation matrix C");
  //   int logPreallocCLoop = 1;  //   PetscLogEventRegister("Prealloc Matrix C", 0, &logPreallocCLoop);
  //   PetscLogEventBegin(logPreallocCLoop, 0, 0, 0, 0);
  //   PetscInt nnz[n]; // Number of non-zeros per row

  //   if (utils::MasterSlave::_rank <= 0)
  //     for (int i = 0; i < polyparams; i++)
  //       nnz[i] = n; // The first rows are dense, except the part in upper right corner

  //   for (const mesh::Vertex& inVertex : inMesh->vertices()) {
  //     if (not inVertex.isOwner())
  //       continue;

  //     int row = inVertex.getGlobalIndex() + polyparams - _matrixC.ownerRange().first; // getGlobalIndex ist absolut, row[] ist zero based
  //     // preciceDebug("Row = " << row << " n = " << n);
  //     nnz[row] = 0;
  
  //     for (mesh::Vertex& vj : inMesh->vertices()) {
  //       distance = inVertex.getCoords() - vj.getCoords();
  //       for (int d = 0; d < dimensions; d++) {
  //         if (_deadAxis[d]) {
  //           distance[d] = 0;
  //         }
  //       }
  //       double coeff = _basisFunction.evaluate(norm2(distance));
  //       if (not tarch::la::equals(coeff, 0.0))
  //         nnz[row]++;
  //     }
  //     // preciceDebug("Reserved for row = " << row << ", reserved = " << nnz[row]);
  //   }
  //   PetscLogEventEnd(logPreallocCLoop, 0, 0, 0, 0);
  //   // -- END PREALLOC LOOP FOR MATRIX C --

  //   ierr = MatMPISBAIJSetPreallocation(_matrixC, 1, 0, nnz, 0, nnz); CHKERRV(ierr);    
  // }
  // -- BEGIN FILL LOOP FOR MATRIX C --
  int logCLoop = 2;
  PetscLogEventRegister("Filling Matrix C", 0, &logCLoop);
  PetscLogEventBegin(logCLoop, 0, 0, 0, 0);
  precice::utils::Event eFillC("Filling Matrix C");
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
    ierr = MatSetValuesLocal(_matrixC, 1, &polyRow, 1, &polyCol, &y, INSERT_VALUES); CHKERRV(ierr);
    int actualDim = 0;
    for (int dim = 0; dim < dimensions; dim++) {
      if (not _deadAxis[dim]) {
        y = inVertex.getCoords()[dim];
        polyRow = 1 + actualDim;
        actualDim++;
        ierr = MatSetValuesLocal(_matrixC, 1, &polyRow, 1, &polyCol, &y, INSERT_VALUES); CHKERRV(ierr);
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
    }

    ierr = MatSetValuesLocal(_matrixC, 1, &row, colNum, colIdx, colVals, INSERT_VALUES); CHKERRV(ierr);
  }
  preciceDebug("Finished filling Matrix C");
  eFillC.stop();
  PetscLogEventEnd(logCLoop, 0, 0, 0, 0);
  // -- END FILL LOOP FOR MATRIX C --

  // Petsc requires that all diagonal entries are set, even if set to zero.
  _matrixC.assemble(MAT_FLUSH_ASSEMBLY);
  petsc::Vector zeros(_matrixC);
  MatDiagonalSet(_matrixC, zeros, ADD_VALUES);

  // Begin assembly here, all assembly is ended at the end of this function.
  ierr = MatAssemblyBegin(_matrixC, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
  
  // -- BEGIN PREALLOC LOOP FOR MATRIX A --
  // {
  //   int localDiagColBegin = _matrixA.ownerRangeColumn().first;
  //   int localDiagColEnd = _matrixA.ownerRangeColumn().second;
  //   preciceDebug("Local Submatrix Rows = " << ownerRangeABegin << " / " << ownerRangeAEnd <<
  //                ", Local Submatrix Cols = " << localDiagColBegin << " / " << localDiagColEnd);

  //   preciceDebug("Begin preallocation matrix A.");
  //   int logPreallocALoop = 3;
  //   PetscLogEventRegister("Prealloc Matrix A", 0, &logPreallocALoop);
  //   PetscLogEventBegin(logPreallocALoop, 0, 0, 0, 0);
  //   PetscInt nnzA[outputSize]; // Number of non-zero entries, off-diagonal
  //   PetscInt DnnzA[outputSize]; // Number of non-zero entries, diagonal
  //   for (int i = 0; i < ownerRangeAEnd - ownerRangeABegin; i++) {
  //     nnzA[i] = 0; DnnzA[i] = 0;
  //     int polyCol = 1;
  //     for (int dim = 0; dim < dimensions; dim++) {
  //       if (not _deadAxis[dim]) {
  //         incPrealloc(nnzA, DnnzA, polyCol, localDiagColBegin, localDiagColEnd);
  //         ++polyCol;
  //       }
  //     }
  //     const mesh::Vertex& oVertex = outMesh->vertices()[i];
  //     for (const mesh::Vertex& inVertex : inMesh->vertices()) {
  //       distance = oVertex.getCoords() - inVertex.getCoords();
  //       for (int d = 0; d < dimensions; d++) {
  //         if (_deadAxis[d])
  //           distance[d] = 0;
  //       }
  //       double coeff = _basisFunction.evaluate(norm2(distance));
  //       if (not tarch::la::equals(coeff, 0.0)) {
  //         const int colPosition = inVertex.getGlobalIndex() + polyparams;
  //         incPrealloc(nnzA, DnnzA, colPosition, localDiagColBegin, localDiagColEnd);
  //       }
  //       preciceDebug("Preallocating     diagonal row " << i << " with " << DnnzA[i] << " elements.");
  //       preciceDebug("Preallocating off-diagonal row " << i << " with " << nnzA[i] << " elements.");
  //     }
  //   }
  //   // -- END PREALLOC LOOP FOR MATRIX A --

  //   // Preallocation: Number of diagonal non-zero entries equals outputSize, since diagonal is full
  //   MatSetOption(_matrixA, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  //   ierr = MatMPIAIJSetPreallocation(_matrixA, 0, DnnzA, 0, nnzA); CHKERRV(ierr);
  // }
  
  // -- BEGIN FILL LOOP FOR MATRIX A --
  preciceDebug("Begin filling matrix A.");
  int logALoop = 4;
  PetscLogEventRegister("Filling Matrix A", 0, &logALoop);
  PetscLogEventBegin(logALoop, 0, 0, 0, 0);
  precice::utils::Event eFillA("Filling Matrix A");

  for (int it = ownerRangeABegin; it < ownerRangeAEnd; it++) {
    // hier colIdx, colVals über inMesh->vertices.count dimensionieren?
    PetscInt colNum = 0;
    PetscInt colIdx[_matrixC.getSize().second];     // holds the columns indices of the entries MatrixC??
    PetscScalar colVals[_matrixC.getSize().second]; // holds the values of the entries
    const mesh::Vertex& oVertex = outMesh->vertices()[it - _matrixA.ownerRange().first];

    // -- SET THE POLYNOM PART OF THE MATRIX --
    PetscInt polyRow = it, polyCol = 0;
    PetscScalar y = 1;
    ierr = MatSetValuesLocal(_matrixA, 1, &polyRow, 1, &polyCol, &y, INSERT_VALUES); CHKERRV(ierr);
    for (int dim = 0; dim < dimensions; dim++) {
      if (not _deadAxis[dim]) {
        y = oVertex.getCoords()[dim];
        polyCol++;
        // preciceDebug("Filling A with polyparams: polyRow = " << polyRow << ", polyCol = " << polyCol << " Preallocation = " << nnzA[polyRow]);
        ierr = MatSetValuesLocal(_matrixA, 1, &polyRow, 1, &polyCol, &y, INSERT_VALUES); CHKERRV(ierr); // das zusammen mit den MatSetValuesLocal unten machen
      }
    }

    // -- SETS THE COEFFICIENTS --
    for (const mesh::Vertex& inVertex : inMesh->vertices()) {
      distance = oVertex.getCoords() - inVertex.getCoords();
      for (int d = 0; d < dimensions; d++) {
        if (_deadAxis[d])
          distance[d] = 0;
      }
      double coeff = _basisFunction.evaluate(norm2(distance));
      if (not tarch::la::equals(coeff, 0.0)) {
        colVals[colNum] = coeff;
        colIdx[colNum] = inVertex.getGlobalIndex() + polyparams;
        colNum++;
      }
    }
    // preciceDebug("Filling A: row = " << it << ", col count = " << colNum);
    ierr = MatSetValuesLocal(_matrixA, 1, &it, colNum, colIdx, colVals, INSERT_VALUES); CHKERRV(ierr);
  }
  preciceDebug("Finished filling Matrix A");
  eFillA.stop();
  PetscLogEventEnd(logALoop, 0, 0, 0, 0);
  // -- END FILL LOOP FOR MATRIX A --

  ierr = MatAssemblyBegin(_matrixA, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
  ierr = MatAssemblyEnd(_matrixC, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
  ierr = MatAssemblyEnd(_matrixA, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
  KSPSetOperators(_solver, _matrixC, _matrixC); CHKERRV(ierr);
  KSPSetTolerances(_solver, _solverRtol, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
  KSPSetInitialGuessNonzero(_solver, PETSC_TRUE); CHKERRV(ierr); // Reuse the results from the last iteration, held in the out vector.
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

  preciceDebug("Number of mallocs for matrix C = " << _matrixC.getInfo(MAT_LOCAL).mallocs);
  preciceDebug("Non-zeros allocated / used / unused for matrix C = " << _matrixC.getInfo(MAT_LOCAL).nz_allocated << " / " << _matrixC.getInfo(MAT_LOCAL).nz_used << " / " << _matrixC.getInfo(MAT_LOCAL).nz_unneeded);
  preciceDebug("Number of mallocs for matrix A = " << _matrixA.getInfo(MAT_LOCAL).mallocs);
  preciceDebug("Non-zeros allocated / used / unused for matrix A = " << _matrixA.getInfo(MAT_LOCAL).nz_allocated << " / " << _matrixA.getInfo(MAT_LOCAL).nz_used << " / " << _matrixA.getInfo(MAT_LOCAL).nz_unneeded);
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
  previousSolution.clear();
  _hasComputedMapping = false;
  // ISmapping destroy??
}

template<typename RADIAL_BASIS_FUNCTION_T>
void PetRadialBasisFctMapping<RADIAL_BASIS_FUNCTION_T>::map(int inputDataID, int outputDataID)
{
  precice::utils::Event e(__func__);
  preciceTrace("map()", inputDataID, outputDataID);

  assertion(_hasComputedMapping);
  assertion(input()->getDimensions() == output()->getDimensions(),
            input()->getDimensions(), output()->getDimensions());

  const std::string constraintName = getConstraint() == CONSERVATIVE ? "conservative" : "consistent";
  preciceInfo(__func__, "Mapping " << input()->data(inputDataID)->getName()
              << " " << constraintName
              << " from " << input()->getName() << " (ID " << input()->getID() << ")"
              << " to " << output()->getName() << " (ID " << output()->getID() << ")");

  PetscErrorCode ierr = 0;
  KSPConvergedReason convReason;
  auto& inValues = input()->data(inputDataID)->values();
  auto& outValues = output()->data(outputDataID)->values();

  int valueDim = input()->data(inputDataID)->getDimensions();
  assertion(valueDim == output()->data(outputDataID)->getDimensions(),
            valueDim, output()->data(outputDataID)->getDimensions());
  int deadDimensions = 0;
  for (int d=0; d<getDimensions(); d++) {
    if (_deadAxis[d]) deadDimensions +=1;
  }
  int polyparams = 1 + getDimensions() - deadDimensions;
  int localPolyparams = utils::MasterSlave::_rank > 0 ? 0 : polyparams; // Set localPolyparams only when root rank

  if (getConstraint() == CONSERVATIVE) {
    petsc::Vector Au(_matrixC, "Au");
    petsc::Vector in(_matrixA, "in");
    
    // Fill input from input data values
    for (int dim=0; dim < valueDim; dim++) {
      preciceDebug("input()->vertices().size() = " << input()->vertices().size());
      for (size_t i = 0; i < input()->vertices().size(); i++ ) {
        int globalIndex = input()->vertices()[i].getGlobalIndex();
        VecSetValue(in, globalIndex, inValues[(i)*valueDim + dim], INSERT_VALUES); // Dies besser als VecSetValuesLocal machen
      }
      in.assemble();

      // Gets the petsc::vector for the given combination of outputData, inputData and dimension
      // If none created yet, create one, based on _matrixC
      petsc::Vector& out = std::get<0>(
        previousSolution.emplace(std::piecewise_construct,
                                 std::forward_as_tuple(inputDataID + outputDataID * 10 + dim * 100),
                                 std::forward_as_tuple(_matrixC, "out"))
        )->second;
      
      ierr = MatMultTranspose(_matrixA, in, Au); CHKERRV(ierr);
      ierr = KSPSolve(_solver, Au, out); CHKERRV(ierr);
      ierr = KSPGetConvergedReason(_solver, &convReason); CHKERRV(ierr);
      if (convReason < 0) {
        KSPView(_solver, PETSC_VIEWER_STDOUT_WORLD);
        preciceError(__func__, "RBF linear system has not converged.");
      }
      
      // petsc::Vector res(_matrixC, "Residual");
      // PetscReal resNorm;
      // MatResidual(_matrixC, Au, out, res);
      // VecNorm(res, NORM_2, &resNorm);
      // res.view();
      // preciceDebug("Residual norm = " << resNorm);
      
      VecChop(out, 1e-9);

      // Copy mapped data to output data values
      const PetscScalar *outArray;
      VecGetArrayRead(out, &outArray);

      int count = 0, ownerCount = 0;
      for (const mesh::Vertex& vertex : output()->vertices()) {
        if (vertex.isOwner()) {
          outValues[count*valueDim + dim] = outArray[ownerCount+localPolyparams];
          ownerCount++;
        }
        else {
          assertion(outValues[count*valueDim + dim] == 0.0);
        }
        count++;
      }
      VecRestoreArrayRead(out, &outArray);
    }
  }
  else { // Map CONSISTENT
    petsc::Vector out(_matrixA, "out");
    petsc::Vector in(_matrixC, "in");
        
    ierr = VecSetLocalToGlobalMapping(in, _ISmapping); CHKERRV(ierr);
    const PetscScalar *vecArray;

    // For every data dimension, perform mapping
    for (int dim=0; dim < valueDim; dim++) {
      // Fill input from input data values
      int count = 0;
      for (const auto& vertex : input()->vertices()) {
        VecSetValueLocal(in, vertex.getGlobalIndex()+polyparams, inValues[count*valueDim + dim], INSERT_VALUES); // evtl. besser als VecSetValuesLocal
        count++;
      }
      in.assemble();
      
      petsc::Vector& p = std::get<0>(  // Save and reuse the solution from the previous iteration
        previousSolution.emplace(std::piecewise_construct,
                                 std::forward_as_tuple(inputDataID + outputDataID * 10 + dim * 100),
                                 std::forward_as_tuple(_matrixC, "p"))
        )->second;
      
      ierr = KSPSolve(_solver, in, p); CHKERRV(ierr);
      ierr = KSPGetConvergedReason(_solver, &convReason); CHKERRV(ierr);
      if (convReason < 0) {
        KSPView(_solver, PETSC_VIEWER_STDOUT_WORLD);
        preciceError(__func__, "RBF linear system has not converged.");
      }
      ierr = MatMult(_matrixA, p, out); CHKERRV(ierr);
      VecChop(out, 1e-9);

      // Copy mapped data to output data values
      ierr = VecGetArrayRead(out, &vecArray);
      int size = out.getLocalSize();
      for (int i=0; i < size; i++) {
        outValues[i*valueDim + dim] = vecArray[i];
      }
      VecRestoreArrayRead(out, &vecArray);
    }
  }
}


template<typename RADIAL_BASIS_FUNCTION_T>
bool PetRadialBasisFctMapping<RADIAL_BASIS_FUNCTION_T>::doesVertexContribute(int vertexID) const
{
  // FIXME: Use a sane calculation here
  // preciceTrace(__func__);

  if (not _basisFunction.hasCompactSupport())
    return true;

  return true;
}

template<typename RADIAL_BASIS_FUNCTION_T>
void PetRadialBasisFctMapping<RADIAL_BASIS_FUNCTION_T>::incPrealloc(PetscInt* diag, PetscInt* offDiag, int pos, int begin, int end)
{
  // Do some optimizatin here: inline, const, ...
  if ((pos < begin) or (pos > end))
    (*offDiag)++; // vertex is off-diagonal
  else
    (*diag)++; // vertex is diagonal
}

}} // namespace precice, mapping

#endif
