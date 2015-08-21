#pragma once
#ifndef PRECICE_NO_PETSC

#include "mapping/Mapping.hpp"
#include "RadialBasisFctMapping.hpp"
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
   * @param constraint [IN] Specifies mapping to be consistent or conservative.
   * @param function [IN] Radial basis function used for mapping.
   * @param solverRtol [IN] Relative tolerance for the linear solver.
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
  virtual bool hasComputedMapping();

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

  // FIXME: Hack to get global index also when MasterSlave mode is not enabled.
  void addGlobalIndex(mesh::PtrMesh &mesh);

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
    preciceCheck(not (xDead && yDead), "setDeadAxis()", "You cannot  "
                 << " choose all axis to be dead for a RBF mapping");
    preciceCheck(not zDead, "setDeadAxis()", "You cannot  "
                 << " dead out the z axis if dimension is set to 2");
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

//  preciceCheck(not utils::MasterSlave::_slaveMode && not utils::MasterSlave::_masterMode,
//               "computeMapping()", "RBF mapping "
//               << "is not yet supported for a participant in master mode");

  using namespace tarch::la;
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

  // Add global indizes to vertices. Benjamin will fix and implement this in a more generic fashion soon
  if (not (utils::MasterSlave::_slaveMode or utils::MasterSlave::_masterMode)) {
    addGlobalIndex(inMesh);
    addGlobalIndex(outMesh);
  }

  // Indizes for Q^T, holding the polynom
  std::vector<int> myIndizes; // vector::reserve for efficiency?
  if (utils::MasterSlave::_rank <= 0)
    for (size_t i = 0; i < polyparams; i++)
      myIndizes.push_back(i);

  // Indizes for the vertices with polyparams offset
  for (const mesh::Vertex& v : inMesh->vertices())
    if (v.isOwner())
      myIndizes.push_back(v.getGlobalIndex() + polyparams);

  auto inputSize = myIndizes.size();
  auto n = inputSize; // polyparams, if on rank 0, are included here

  if (utils::MasterSlave::_rank <= 0)
    inputSize -= polyparams;

  auto outputSize = outMesh->vertices().size();

  assertion1(inputSize >= 1 + polyparams, inputSize);

  PetscErrorCode ierr = 0;

  _matrixC.reset();
  ierr = MatSetType(_matrixC.matrix, MATSBAIJ); CHKERRV(ierr); // create symmetric, block sparse matrix.
  preciceDebug("Set matrix C to size " << n);
  ierr = MatSetSizes(_matrixC.matrix, n, n, PETSC_DECIDE, PETSC_DECIDE); CHKERRV(ierr);
  // ierr = MatSetOption(_matrixC.matrix, MAT_SYMMETRY_ETERNAL, PETSC_TRUE); CHKERRV(ierr);
  ierr = MatSetUp(_matrixC.matrix); CHKERRV(ierr);

  _matrixA.reset();
  ierr = MatSetType(_matrixA.matrix, MATAIJ); CHKERRV(ierr); // create sparse matrix.
  ierr = MatSetSizes(_matrixA.matrix, outputSize, n, PETSC_DECIDE, PETSC_DECIDE); CHKERRV(ierr);
  ierr = MatSetUp(_matrixA.matrix); CHKERRV(ierr);

  PetscInt range_start, range_end;
  MatGetOwnershipRange(_matrixA.matrix, &range_start, &range_end);

  KSPReset(_solver);

  IS ISlocal, ISglobal, ISidentity, ISidentityGlobal;
  ISLocalToGlobalMapping ISidentityMapping;
  ierr = ISCreateGeneral(PETSC_COMM_WORLD, myIndizes.size(), myIndizes.data(), PETSC_COPY_VALUES, &ISlocal); CHKERRV(ierr);
  ierr = ISAllGather(ISlocal, &ISglobal); CHKERRV(ierr);
  ierr = ISLocalToGlobalMappingCreateIS(ISglobal, &_ISmapping); CHKERRV(ierr);
  ierr = MatSetLocalToGlobalMapping(_matrixC.matrix, _ISmapping, _ISmapping); CHKERRV(ierr);
  // Create an identy mapping and use that for the columns of matrixA.
  ierr = ISCreateStride(PETSC_COMM_WORLD, _matrixA.ownerRange().second - _matrixA.ownerRange().first, _matrixA.ownerRange().first, 1, &ISidentity); CHKERRV(ierr);
  ierr = ISAllGather(ISidentity, &ISidentityGlobal); CHKERRV(ierr);
  ierr = ISLocalToGlobalMappingCreateIS(ISidentityGlobal, &ISidentityMapping); CHKERRV(ierr);
  ierr = MatSetLocalToGlobalMapping(_matrixA.matrix, ISidentityMapping, _ISmapping); CHKERRV(ierr);

  int myRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  // if (myRank == 0) {
  // cout << "== ISidentityGlobal ==" << endl;
  // ISView(ISidentityGlobal, PETSC_VIEWER_STDOUT_SELF);
  // ISView(ISglobal, PETSC_VIEWER_STDOUT_SELF);
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
  PetscInt colIdx[n];     // holds the columns indices of the entries
  PetscScalar colVals[n]; // holds the values of the entries
  for (const mesh::Vertex& iVertex : inMesh->vertices()) {
    PetscInt polyRow = 0, polyCol = polyparams+iVertex.getGlobalIndex();
    PetscScalar y = 1;
    ierr = MatSetValuesLocal(_matrixC.matrix, 1, &polyRow, 1, &polyCol, &y, INSERT_VALUES); CHKERRV(ierr);
    int actualDim = 0;
    for (int dim = 0; dim < dimensions; dim++) {
      if (not _deadAxis[dim]) {
        y = iVertex.getCoords()[dim];
        polyRow = 1 + actualDim;
        actualDim++;
        ierr = MatSetValuesLocal(_matrixC.matrix, 1, &polyRow, 1, &polyCol, &y, INSERT_VALUES); CHKERRV(ierr);
      }
    }

    PetscInt colNum = 0;  // holds the number of columns
    for (size_t j=iVertex.getID(); j < inputSize; j++) {
      distance = iVertex.getCoords() - inMesh->vertices()[j].getCoords();
      for (int d = 0; d < dimensions; d++) {
        if (_deadAxis[d]) {
          distance[d] = 0;
        }
      }
      double coeff = _basisFunction.evaluate(norm2(distance));
      if (not equals(coeff, 0.0)) {
        colVals[colNum] = coeff;
        colIdx[colNum] = inMesh->vertices()[j].getGlobalIndex() + polyparams; // column of entry is the globaIndex
        colNum++;
      }
      #ifdef Asserts
      if (coeff == std::numeric_limits<double>::infinity()) {
        preciceError("computeMapping()", "C matrix element has value inf. "
                     << "i = " << i << ", j = " << j
                     << ", coords i = " << iVertex.getCoords() << ", coords j = "
                     << inMesh->vertices()[j].getCoords() << ", dist = "
                     << distance << ", norm2 = " << norm2(distance) << ", rbf = "
                     << coeff
                     << ", rbf type = " << typeid(_basisFunction).name());
      }
      # endif
    }
    // colVals[colNum] = 1;
    // colIdx[colNum] = inputSize;
    // colNum++;
    // int actualDim = 0;
    // for (int dim=0; dim < dimensions; dim++) {
    //   if (not _deadAxis[dim]) {
    //     colVals[colNum] = iVertex.getCoords()[dim];
    //     colIdx[colNum] = inputSize+1+actualDim;
    //     colNum++;
    //     actualDim++;
    //   }
    // }
    int row = iVertex.getGlobalIndex() + polyparams;
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
  // ierr = MatAssemblyBegin(_matrixC.matrix, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);

  _matrixC.assemble();

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
  for (const mesh::Vertex& oVertex : outMesh->vertices()) {
    // _matrixA.assemble();
    // _matrixA.view();
    PetscInt polyRow = oVertex.getGlobalIndex(), polyCol = 0;
    PetscScalar y = 1;
    ierr = MatSetValuesLocal(_matrixA.matrix, 1, &polyRow, 1, &polyCol, &y, INSERT_VALUES); CHKERRV(ierr);
    for (int dim = 0; dim < dimensions; dim++) {
      if (not _deadAxis[dim]) {
        y = oVertex.getCoords()[dim];
        polyCol++;
        ierr = MatSetValuesLocal(_matrixA.matrix, 1, &polyRow, 1, &polyCol, &y, INSERT_VALUES); CHKERRV(ierr);
      }
    }

    PetscInt colNum = 0;
    int j = 0;
    for (const mesh::Vertex& jVertex : inMesh->vertices()) {
      distance = oVertex.getCoords() - jVertex.getCoords();
      for (int d = 0; d < dimensions; d++) {
        if (_deadAxis[d])
          distance[d] = 0;
      }
      double coeff = _basisFunction.evaluate(norm2(distance));
      if ( not equals(coeff, 0.0)) {
        colVals[colNum] = coeff;
        colIdx[colNum] = jVertex.getGlobalIndex() + polyparams;
        colNum++;
      }
      #     ifdef Asserts
      if (coeff == std::numeric_limits<double>::infinity()){
        preciceError("computeMapping()", "A matrix element has value inf. "
                     << "i = " << i << ", j = " << j
                     << ", coords i = " << oVertex.getCoords() << ", coords j = "
                     << jVertex.getCoords() << ", dist = "
                     << distance << ", norm2 = " << norm2(distance) << ", rbf = "
                     << coeff
                     << ", rbf type = " << typeid(_basisFunction).name());
      }
      #     endif
      j++;
    }
    // ierr = MatSetValue(_matrixA.matrix, i, inputSize, 1.0, INSERT_VALUES); CHKERRV(ierr);
    // colVals[colNum] = 1;
    // colIdx[colNum] = inputSize;
    // colNum++;
    // int actualDim = 0;
    // for (int dim=0; dim < dimensions; dim++) {
    //   if (not _deadAxis[dim]) {
    //     colVals[colNum] = iVertex.getCoords()[dim];
    //     colIdx[colNum] = inputSize+1+actualDim;
    //     colNum++;
    //     actualDim++;
    //   }
    // }
    ierr = MatSetValuesLocal(_matrixA.matrix, 1, &i, colNum, colIdx, colVals, INSERT_VALUES); CHKERRV(ierr);
    i++;
  }
  PetscLogEventEnd(logALoop, 0, 0, 0, 0);
  // -- END FILL LOOP FOR MATRIX A --

  ierr = MatAssemblyBegin(_matrixA.matrix, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
  // ierr = MatAssemblyEnd(_matrixC.matrix, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
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
  // _matrixC.view();
  _matrixA.view();
}

template<typename RADIAL_BASIS_FUNCTION_T>
bool PetRadialBasisFctMapping<RADIAL_BASIS_FUNCTION_T>:: hasComputedMapping()
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
      int size = in.getSize();
      // for (int i=0; i < size; i++) { // Fill input data values
      // in.setValue(i, inValues[i*valueDim + dim]);
      // }
      for (int i=0; i < size; i++) {
        int globalIndex = input()->vertices()[i].getGlobalIndex();
        // Dies besser als VecSetValuesLocal machen
        VecSetValueLocal(in.vector, globalIndex, inValues[i*valueDim + dim], INSERT_VALUES);
      }
      in.assemble();

      ierr = MatMultTranspose(_matrixA.matrix, in.vector, Au.vector); CHKERRV(ierr);
      ierr = KSPSolve(_solver, Au.vector, out.vector); CHKERRV(ierr);
      ierr = KSPGetConvergedReason(_solver, &convReason); CHKERRV(ierr);
      if (convReason < 0) {
        preciceError(__func__, "RBF linear system has not converged.");
      }
      VecChop(out.vector, 1e-9);
      // Copy mapped data to output data values
      PetscScalar *outArray;
      ierr = VecGetArray(out.vector, &outArray);
      size = out.getSize();
      for (int i=0; i < size-polyparams; i++){
        outValues[i*valueDim + dim] = outArray[i+polyparams]; // hier noch das index set beachten?
      }
      VecRestoreArray(out.vector, &outArray);
    }
    mappingIndex++;
  }
  else { // Map consistent
    preciceDebug("Map consistent");
    petsc::Vector p(_matrixC, "p");
    petsc::Vector in(_matrixC, "in");
    ierr = VecSetLocalToGlobalMapping(in.vector, _ISmapping); CHKERRV(ierr);
    petsc::Vector out(_matrixA, "out");
    PetscScalar *vecArray;

    // For every data dimension, perform mapping
    for (int dim=0; dim < valueDim; dim++){
      // Fill input from input data values (last polyparams entries remain zero)
      // ierr = VecGetArray(in.vector, &vecArray);
      int size  = in.getSize();
      for (int i=polyparams; i < size; i++){
        int globalIndex = input()->vertices()[i-polyparams].getGlobalIndex();
        // vecArray[globalIndex] = inValues[(i-polyparams)*valueDim + dim];
        // Dies besser als VecSetValuesLocal machen
        VecSetValueLocal(in.vector, globalIndex+polyparams, inValues[(i-polyparams)*valueDim + dim], INSERT_VALUES);
      }
      // VecRestoreArray(in.vector, &vecArray);
      in.assemble();

      ierr = KSPSolve(_solver, in.vector, p.vector); CHKERRV(ierr);
      ierr = KSPGetConvergedReason(_solver, &convReason); CHKERRV(ierr);
      if (convReason < 0) {
        preciceError(__func__, "RBF linear system has not converged.");
      }
      ierr = MatMult(_matrixA.matrix, p.vector, out.vector); CHKERRV(ierr);
      VecChop(out.vector, 1e-9);
      // Copy mapped data to output data values
      ierr = VecGetArray(out.vector, &vecArray);
      size = out.getSize();
      for (int i=0; i < size; i++) {
        outValues[i*valueDim + dim] = vecArray[i];
      }
      VecRestoreArray(out.vector, &vecArray);
    }
  }
}

template<typename RADIAL_BASIS_FUNCTION_T>
void PetRadialBasisFctMapping<RADIAL_BASIS_FUNCTION_T>::addGlobalIndex(mesh::PtrMesh &mesh)
{
  preciceTrace(__func__);
  preciceDebug("Manually adding indizes to vertices.");
  size_t i = 0;
  for (mesh::Vertex& v : mesh->vertices()) {
    v.setGlobalIndex(i);
    // v.setGlobalIndex(v.getID());
    i++;
  }
}


}} // namespace precice, mapping

#endif
