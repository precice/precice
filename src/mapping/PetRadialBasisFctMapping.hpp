#pragma once
#ifndef PRECICE_NO_PETSC

#include "mapping/Mapping.hpp"
#include "RadialBasisFctMapping.hpp"
#include "tarch/la/DynamicVector.h"
#include "utils/MasterSlave.hpp"
#include <limits>
#include <typeinfo>

#include "petnum.hpp"
#include "petscmat.h"
#include "petscksp.h"
#include "petsclog.h"

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

  double _solverRtol;

  /// true if the mapping along some axis should be ignored
  bool* _deadAxis;

  /// Deletes all dead directions from fullVector and returns a vector of reduced dimensionality.
  // utils::DynVector reduceVector(const utils::DynVector& fullVector);

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
  _matrixC(PETSC_COMM_SELF, "C"),
  _matrixA(PETSC_COMM_SELF, "A"),
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
    
  KSPCreate(PETSC_COMM_SELF, &_solver);
}

template<typename RADIAL_BASIS_FUNCTION_T>
PetRadialBasisFctMapping<RADIAL_BASIS_FUNCTION_T>::~PetRadialBasisFctMapping()
{
  delete[] _deadAxis;
  KSPDestroy(&_solver);
}



template<typename RADIAL_BASIS_FUNCTION_T>
void PetRadialBasisFctMapping<RADIAL_BASIS_FUNCTION_T>::computeMapping()
{
  preciceTrace("computeMapping()");

  preciceCheck(not utils::MasterSlave::_slaveMode && not utils::MasterSlave::_masterMode,
               "computeMapping()", "RBF mapping "
               << "is not yet supported for a participant in master mode");

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
  int inputSize = (int)inMesh->vertices().size();
  int outputSize = (int)outMesh->vertices().size();
  int deadDimensions = 0;
  for (int d=0; d<dimensions; d++) {
    if (_deadAxis[d]) deadDimensions +=1;
  }

  int polyparams = 1 + dimensions - deadDimensions;
  PetscErrorCode ierr = 0;
  assertion1(inputSize >= 1 + polyparams, inputSize);
  int n = inputSize + polyparams; // Add linear polynom degrees

  _matrixC.reset(); 
  ierr = MatSetType(_matrixC.matrix, MATSBAIJ); CHKERRV(ierr); // create symmetric, block sparse matrix.
  ierr = MatSetSizes(_matrixC.matrix, PETSC_DECIDE, PETSC_DECIDE, n, n); CHKERRV(ierr);
  ierr = MatSetOption(_matrixC.matrix, MAT_SYMMETRY_ETERNAL, PETSC_TRUE); CHKERRV(ierr);

  _matrixA.reset();
  ierr = MatSetType(_matrixA.matrix, MATAIJ); CHKERRV(ierr); // create sparse matrix.
  ierr = MatSetSizes(_matrixA.matrix, PETSC_DECIDE, PETSC_DECIDE, outputSize, n); CHKERRV(ierr);

  KSPReset(_solver);

  int i = 0;
  utils::DynVector distance(dimensions);

  // We do preallocating of the matrices C and A. That means we traverse the input data once, just
  // to know where we have entries in the sparse matrix. This information petsc can use to
  // preallocate the matrix. In the second phase we actually fill the matrix.

  // -- BEGIN PREALLOC LOOP FOR MATRIX C --
  int logPreallocCLoop = 1;
  PetscLogEventRegister("Prealloc Matrix C", 0, &logPreallocCLoop);
  PetscLogEventBegin(logPreallocCLoop, 0, 0, 0, 0);
  PetscInt nnz[n]; // Number of non-zeros per row
  unsigned int totalNNZ = 0; // Total number of non-zeros in matrix
  for (const mesh::Vertex& iVertex : inMesh->vertices()) {
    nnz[i] = 0;
    for (int j=iVertex.getID(); j < inputSize; j++) {
      if (i == j) {
        nnz[i]++; // Since we need to set at least zeros on the main diagonal, we have an entry
                  // there. No further test necessary.
        continue;
      }
      distance = iVertex.getCoords() - inMesh->vertices()[j].getCoords();
      double coeff = _basisFunction.evaluate(norm2(distance));
      if ( not equals(coeff, 0.0)) {
        nnz[i]++;
      }
    }
    nnz[i] += dimensions + 1 - deadDimensions; // coefficients of the linear polynom
    totalNNZ += nnz[i];
    i++;
  }
  for (int r = inputSize; r < n; r++) {
    nnz[r] = 1; // iterate through the lower part, no coefficients there, but 0 on the main diagonal
  }
  PetscLogEventEnd(logPreallocCLoop, 0, 0, 0, 0);
  i = 0;
  // -- END PREALLOC LOOP FOR MATRIX C --
  
  ierr = MatSeqSBAIJSetPreallocation(_matrixC.matrix, 1, PETSC_DEFAULT, nnz);
  
  // -- BEGIN FILL LOOP FOR MATRIX C --
  int logCLoop = 2;
  PetscLogEventRegister("Filling Matrix C", 0, &logCLoop);
  PetscLogEventBegin(logCLoop, 0, 0, 0, 0);
  // We collect entries for each row and set them blockwise using MatSetValues.
  PetscInt colIdx[n];     // holds the columns indices of the entries
  PetscScalar colVals[n]; // holds the values of the entries
  for (const mesh::Vertex& iVertex : inMesh->vertices()) {
    PetscInt colNum = 0;  // holds the number of columns
    for (int j=iVertex.getID(); j < inputSize; j++) {
      distance = iVertex.getCoords() - inMesh->vertices()[j].getCoords();
      for (int d = 0; d < dimensions; d++) {
        if (_deadAxis[d]) {
          distance[d] = 0;
        }
      }
      double coeff = _basisFunction.evaluate(norm2(distance));
      if ( not equals(coeff, 0.0)) {
        colVals[colNum] = coeff;
        colIdx[colNum] = j;
        colNum++;
      }
#     ifdef Asserts
      if (coeff == std::numeric_limits<double>::infinity()) {
        preciceError("computeMapping()", "C matrix element has value inf. "
                     << "i = " << i << ", j = " << j
                     << ", coords i = " << iVertex.getCoords() << ", coords j = "
                     << inMesh->vertices()[j].getCoords() << ", dist = "
                     << distance << ", norm2 = " << norm2(distance) << ", rbf = "
                     << coeff
                     << ", rbf type = " << typeid(_basisFunction).name());
      }
#     endif
    }
    colVals[colNum] = 1;
    colIdx[colNum] = inputSize;
    colNum++;
    int actualDim = 0;
    for (int dim=0; dim < dimensions; dim++) {
      if (not _deadAxis[dim]) {
        colVals[colNum] = iVertex.getCoords()[dim];
        colIdx[colNum] = inputSize+1+actualDim;
        colNum++;
        actualDim++;
      }
    }
    ierr = MatSetValues(_matrixC.matrix, 1, &i, colNum, colIdx, colVals, INSERT_VALUES); CHKERRV(ierr);
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
  int logPreallocALoop = 3;
  PetscLogEventRegister("Prealloc Matrix A", 0, &logPreallocALoop);
  PetscLogEventBegin(logPreallocALoop, 0, 0, 0, 0);
  PetscInt nnzA[outputSize];
  i = 0;
  for (const mesh::Vertex& iVertex : outMesh->vertices()) {
    nnzA[i] = 0;
    for (const mesh::Vertex& jVertex : inMesh->vertices()) {
      distance = iVertex.getCoords() - jVertex.getCoords();
      double coeff = _basisFunction.evaluate(norm2(distance));
      if ( not equals(coeff, 0.0)) {
        nnzA[i]++;
      }
    }
    nnzA[i] += dimensions + 1 - deadDimensions;
    i++;
  }
  PetscLogEventEnd(logPreallocALoop, 0, 0, 0, 0);
  // -- END PREALLOC LOOP FOR MATRIX A --

  ierr = MatSeqAIJSetPreallocation(_matrixA.matrix, PETSC_DEFAULT, nnzA); CHKERRV(ierr);

  // -- BEGIN FILL LOOP FOR MATRIX A --
  int logALoop = 4;
  PetscLogEventRegister("Filling Matrix A", 0, &logALoop);
  PetscLogEventBegin(logALoop, 0, 0, 0, 0);
  i = 0;
  for (const mesh::Vertex& iVertex : outMesh->vertices()) {
    PetscInt colNum = 0;
    int j = 0;
    for (const mesh::Vertex& jVertex : inMesh->vertices()) {
      distance = iVertex.getCoords() - jVertex.getCoords();
      for (int d = 0; d < dimensions; d++) {
        if (_deadAxis[d])
          distance[d] = 0;
      }
      double coeff = _basisFunction.evaluate(norm2(distance));
      if ( not equals(coeff, 0.0)) {
        colVals[colNum] = coeff;
        colIdx[colNum] = j;
        colNum++;
      }
#     ifdef Asserts
      if (coeff == std::numeric_limits<double>::infinity()){
        preciceError("computeMapping()", "A matrix element has value inf. "
                     << "i = " << i << ", j = " << j
                     << ", coords i = " << iVertex.getCoords() << ", coords j = "
                     << jVertex.getCoords() << ", dist = "
                     << distance << ", norm2 = " << norm2(distance) << ", rbf = "
                     << coeff
                     << ", rbf type = " << typeid(_basisFunction).name());
      }
#     endif
      j++;
    }
    ierr = MatSetValue(_matrixA.matrix, i, inputSize, 1.0, INSERT_VALUES); CHKERRV(ierr);
    int actualDim = 0;
    for (int dim=0; dim < dimensions; dim++) {
      if (not _deadAxis[dim]) {
        colVals[colNum] = iVertex.getCoords()[dim];
        colIdx[colNum] = inputSize+1+actualDim;
        colNum++;
        actualDim++;
      }
    }
    ierr = MatSetValues(_matrixA.matrix, 1, &i, colNum, colIdx, colVals, INSERT_VALUES); CHKERRV(ierr);
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

  if (totalNNZ > 20*n) {
    preciceDebug("Using Cholesky decomposition as direct solver for dense matrix.");
    PC prec;
    KSPSetType(_solver, KSPPREONLY);
    KSPGetPC(_solver, &prec);
    PCSetType(prec, PCCHOLESKY);
    PCFactorSetShiftType(prec, MAT_SHIFT_NONZERO);
  }      
  
  _hasComputedMapping = true;
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

    for (int dim=0; dim < valueDim; dim++) {
      int size = in.getSize();
      for (int i=0; i < size; i++) { // Fill input data values
        in.setValue(i, inValues[i*valueDim + dim]);
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
        outValues[i*valueDim + dim] = outArray[i];
      }
      VecRestoreArray(out.vector, &outArray);
    }
    mappingIndex++;
  }
  else { // Map consistent
    preciceDebug("Map consistent");
    petsc::Vector p(_matrixC, "p");
    petsc::Vector in(_matrixC, "in");
    petsc::Vector out(_matrixA, "out");
    PetscScalar *vecArray;

    // For every data dimension, perform mapping
    for (int dim=0; dim < valueDim; dim++){
      // Fill input from input data values (last polyparams entries remain zero)
      ierr = VecGetArray(in.vector, &vecArray);
      int size  = in.getSize();
      for (int i=0; i < size - polyparams; i++){
        vecArray[i] = inValues[i*valueDim + dim];
      }
      VecRestoreArray(in.vector, &vecArray);
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

}} // namespace precice, mapping

#endif
