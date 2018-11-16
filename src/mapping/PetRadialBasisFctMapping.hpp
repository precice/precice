#pragma once
#ifndef PRECICE_NO_PETSC

#include "mapping/Mapping.hpp"

#include <map>
#include <numeric>

#include "versions.hpp"
#include "mesh/RTree.hpp"
#include "math/math.hpp"
#include "impl/BasisFunctions.hpp"
#include "config/MappingConfiguration.hpp"
#include "utils/Petsc.hpp"
namespace petsc = precice::utils::petsc;
#include "utils/EventTimings.hpp"

#include "petscmat.h"
#include "petscksp.h"

// Forward declaration to friend the boost test struct
namespace MappingTests {
namespace PetRadialBasisFunctionMapping {
namespace Serial {
struct SolutionCaching;
}}}

namespace precice {
extern bool testMode;
extern bool syncMode;
}

namespace precice {
namespace mapping {

namespace tests {
class PetRadialBasisFctMappingTest; // Forward declaration to friend the class
}

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
   * @param[in] dimensions Dimensionality of the meshes
   * @param[in] function Radial basis function used for mapping.
   * @param[in] xDead, yDead, zDead Deactivates mapping along an axis
   * @param[in] solverRtol Relative tolerance for the linear solver.
   * @param[in] polynomial Type of polynomial augmentation
   * @param[in] preallocation Sets kind of preallocation of matrices.
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
    double                  solverRtol = 1e-9,
    Polynomial              polynomial = Polynomial::SEPARATE,
    Preallocation           preallocation = Preallocation::TREE);

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

  friend struct MappingTests::PetRadialBasisFunctionMapping::Serial::SolutionCaching;

  virtual void tagMeshFirstRound() override;

  virtual void tagMeshSecondRound() override;


private:
  mutable logging::Logger _log{"mapping::PetRadialBasisFctMapping"};

  bool _hasComputedMapping = false;

  /// Radial basis function type used in interpolation.
  RADIAL_BASIS_FUNCTION_T _basisFunction;

  /// Interpolation system matrix. Evaluated basis function on the input mesh
  petsc::Matrix _matrixC;

  /// Vandermonde Matrix for linear polynomial, constructed from vertices of the input mesh
  petsc::Matrix _matrixQ;

  /// Interpolation evaluation matrix. Evaluated basis function on the output mesh
  petsc::Matrix _matrixA;

  /// Coordinates of the output mesh to evaluate the separated polynomial
  petsc::Matrix _matrixV;

  petsc::Vector rescalingCoeffs;

  /// Used to solve matrixC for the RBF weighting factors
  petsc::KSPSolver _solver;

  /// Used to solve the under-determined system for the separated polynomial.
  petsc::KSPSolver _QRsolver;

  ISLocalToGlobalMapping _ISmapping;

  const double _solverRtol;

  /// true if the mapping along some axis should be ignored
  bool* _deadAxis;

  /// Toggles the use of the additonal polynomial
  Polynomial _polynomial;

  /// Toggles use of rescaled basis functions, only active when Polynomial == SEPARATE
  const bool useRescaling = true;

  /// Number of coefficients for the integrated polynomial. Depends on dimension and number of dead dimensions
  size_t polyparams;

  /// Number of coefficients for the separated polynomial. Depends on dimension and number of dead dimensions
  size_t sepPolyparams;

  /// Equals to polyparams if rank == 0. Is 0 everywhere else
  size_t localPolyparams;

  /// Caches the solution from the previous iteration, used as starting value for current iteration
  std::map<unsigned int, petsc::Vector> previousSolution;

  /// Prints an INFO about the current mapping
  void printMappingInfo(int inputDataID, int dim) const;

  /// Toggles use of preallocation for matrix C and A
  const Preallocation _preallocation;

  void estimatePreallocationMatrixC(int rows, int cols, mesh::PtrMesh mesh);

  void estimatePreallocationMatrixA(int rows, int cols, mesh::PtrMesh mesh);

  void computePreallocationMatrixC(const mesh::PtrMesh inMesh);

  void computePreallocationMatrixA(const mesh::PtrMesh inMesh, const mesh::PtrMesh outMesh);

  std::vector<std::vector<std::pair<int, double>>> savedPreallocationMatrixC(const mesh::PtrMesh inMesh);

  std::vector<std::vector<std::pair<int, double>>> savedPreallocationMatrixA(const mesh::PtrMesh inMesh, const mesh::PtrMesh outMesh);

  std::vector<std::vector<std::pair<int, double>>> bgPreallocationMatrixC(const mesh::PtrMesh inMesh);

  std::vector<std::vector<std::pair<int, double>>> bgPreallocationMatrixA(const mesh::PtrMesh inMesh, const mesh::PtrMesh outMesh);

};

// --------------------------------------------------- HEADER IMPLEMENTATIONS

template<typename RADIAL_BASIS_FUNCTION_T>
PetRadialBasisFctMapping<RADIAL_BASIS_FUNCTION_T>::PetRadialBasisFctMapping
(
  Constraint              constraint,
  int                     dimensions,
  RADIAL_BASIS_FUNCTION_T function,
  bool                    xDead,
  bool                    yDead,
  bool                    zDead,
  double                  solverRtol,
  Polynomial              polynomial,
  Preallocation           preallocation)
  :
  Mapping ( constraint, dimensions ),
  _basisFunction ( function ),
  _matrixC("C"),
  _matrixQ("Q"),
  _matrixA("A"),
  _matrixV("V"),
  _solver("Coefficient Solver"),
  _QRsolver("QR Solver"),
  _ISmapping(nullptr),
  _solverRtol(solverRtol),
  _polynomial(polynomial),
  _preallocation(preallocation)
{
  setInputRequirement(Mapping::MeshRequirement::VERTEX);
  setOutputRequirement(Mapping::MeshRequirement::VERTEX);
  _deadAxis = new bool[dimensions];

  if (getDimensions()==2) {
    _deadAxis[0] = xDead;
    _deadAxis[1] = yDead;
    CHECK(not (xDead && yDead), "You cannot choose all axis to be dead for a RBF mapping");
    CHECK(not zDead, "You cannot dead out the z axis if dimension is set to 2");
  }
  else if (getDimensions()==3) {
    _deadAxis[0] = xDead;
    _deadAxis[1] = yDead;
    _deadAxis[2] = zDead;
    CHECK(not (xDead && yDead && zDead), "You cannot choose all axis to be dead for a RBF mapping");
  }
  else {
    assertion(false);
  }

  // Count number of dead dimensions
  int deadDimensions = 0;
  for (int d = 0; d < dimensions; d++) {
    if (_deadAxis[d]) deadDimensions +=1;
  }
  polyparams =    (_polynomial == Polynomial::ON      ) ? 1 + dimensions - deadDimensions : 0;
  sepPolyparams = (_polynomial == Polynomial::SEPARATE) ? 1 + dimensions - deadDimensions : 0;
  localPolyparams = utils::Parallel::getProcessRank() > 0 ? 0 : polyparams;
}

template<typename RADIAL_BASIS_FUNCTION_T>
PetRadialBasisFctMapping<RADIAL_BASIS_FUNCTION_T>::~PetRadialBasisFctMapping()
{
  delete[] _deadAxis;
  petsc::destroy(&_ISmapping);
}


template<typename RADIAL_BASIS_FUNCTION_T>
void PetRadialBasisFctMapping<RADIAL_BASIS_FUNCTION_T>::computeMapping()
{
  TRACE();
  precice::utils::Event e("map.pet.computeMapping.From" + input()->getName() + "To"+ output()->getName(), precice::syncMode);

  clear();

  if (_polynomial == Polynomial::ON) {
    DEBUG("Using integrated polynomial.");
  }
  if (_polynomial == Polynomial::OFF) {
    DEBUG("Using no polynomial.");
  }
  if (_polynomial == Polynomial::SEPARATE) {
    DEBUG("Using seperated polynomial.");
  }

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

  // Indizes that are used to build the Petsc Index set
  std::vector<PetscInt> myIndizes;

  // Indizes for Q^T, holding the polynomial
  if (utils::Parallel::getProcessRank() <= 0) // Rank 0 or not in MasterSlave mode
    for (size_t i = 0; i < polyparams; i++)
      myIndizes.push_back(i); // polyparams reside in the first rows (which are always on rank 0)

  // Indizes for the vertices with polyparams offset
  for (const mesh::Vertex& v : inMesh->vertices())
    if (v.isOwner())
      myIndizes.push_back(v.getGlobalIndex() + polyparams);

  auto n = myIndizes.size(); // polyparams, if on rank 0, are included here

  auto outputSize = outMesh->vertices().size();

  PetscErrorCode ierr = 0;
  DEBUG("inMesh->vertices().size() = " << inMesh->vertices().size());
  DEBUG("outMesh->vertices().size() = " << outMesh->vertices().size());

  // Matrix C: Symmetric, sparse matrix with n x n local size.
  _matrixC.init(n, n, PETSC_DETERMINE, PETSC_DETERMINE, MATSBAIJ);
  DEBUG("Set matrix C to local size " << n << " x " << n);
  ierr = MatSetOption(_matrixC, MAT_SYMMETRIC, PETSC_TRUE); CHKERRV(ierr);
  ierr = MatSetOption(_matrixC, MAT_SYMMETRY_ETERNAL, PETSC_TRUE); CHKERRV(ierr);

  // Matrix Q: Dense, holds the input mesh for the polynomial if set to SEPERATE. Zero size otherwise
  _matrixQ.init(n, PETSC_DETERMINE, PETSC_DETERMINE, sepPolyparams, MATDENSE);
  DEBUG("Set matrix Q to local size " << n << " x " << sepPolyparams);

  // Matrix V: Dense, holds the output mesh for polynomial if set to SEPERATE. Zero size otherwise
  _matrixV.init(outputSize, PETSC_DETERMINE, PETSC_DETERMINE, sepPolyparams, MATDENSE);
  DEBUG("Set matrix V to local size " << outputSize << " x " << sepPolyparams);

  // Matrix A: Sparse matrix with outputSize x n local size.
  _matrixA.init(outputSize, n, PETSC_DETERMINE, PETSC_DETERMINE, MATAIJ);
  DEBUG("Set matrix A to local size " << outputSize << " x " << n);

  auto const ownerRangeABegin = _matrixA.ownerRange().first;
  auto const ownerRangeAEnd = _matrixA.ownerRange().second;

  IS ISlocal, ISlocalInv, ISglobal, ISidentity, ISidentityGlobal, ISpolyparams;
  ISLocalToGlobalMapping ISidentityMapping, ISpolyparamsMapping;
  // Create an index set which maps myIndizes to continous chunks of matrix rows.
  ierr = ISCreateGeneral(utils::Parallel::getGlobalCommunicator(), myIndizes.size(), myIndizes.data(), PETSC_COPY_VALUES, &ISlocal); CHKERRV(ierr);
  ierr = ISSetPermutation(ISlocal); CHKERRV(ierr);
  ierr = ISInvertPermutation(ISlocal, myIndizes.size(), &ISlocalInv); CHKERRV(ierr);
  ierr = ISAllGather(ISlocalInv, &ISglobal); CHKERRV(ierr); // Gather the IS from all processors
  ierr = ISLocalToGlobalMappingCreateIS(ISglobal, &_ISmapping); CHKERRV(ierr); // Make it a mapping

  // Create an identity mapping and use that for the rows of matrixA.
  ierr = ISCreateStride(utils::Parallel::getGlobalCommunicator(), ownerRangeAEnd - ownerRangeABegin, ownerRangeABegin, 1, &ISidentity); CHKERRV(ierr);
  ierr = ISSetIdentity(ISidentity); CHKERRV(ierr);
  ierr = ISAllGather(ISidentity, &ISidentityGlobal); CHKERRV(ierr);
  ierr = ISLocalToGlobalMappingCreateIS(ISidentityGlobal, &ISidentityMapping); CHKERRV(ierr);

  // Create another identity mapping for the polynomial parameters.
  // That is not needed, but its set, so we can exchange matrixV and matrixA when filling polyparams,
  // i.e. use VecSetValuesLocal on both
  ierr = ISCreateStride(PETSC_COMM_SELF, sepPolyparams, 0, 1, &ISpolyparams); CHKERRV(ierr);
  ierr = ISSetIdentity(ISpolyparams); CHKERRV(ierr);
  ierr = ISLocalToGlobalMappingCreateIS(ISpolyparams, &ISpolyparamsMapping); CHKERRV(ierr);

  // Sets the mappings on the matrices
  ierr = MatSetLocalToGlobalMapping(_matrixC, _ISmapping, _ISmapping); CHKERRV(ierr); // Set mapping for rows and cols
  ierr = MatSetLocalToGlobalMapping(_matrixA, ISidentityMapping, _ISmapping); CHKERRV(ierr); // Set mapping only for cols, use identity for rows
  ierr = MatSetLocalToGlobalMapping(_matrixQ, _ISmapping, ISpolyparamsMapping); CHKERRV(ierr);
  ierr = MatSetLocalToGlobalMapping(_matrixV, ISidentityMapping, ISpolyparamsMapping); CHKERRV(ierr);

  // Destroy all local index sets and mappings
  ierr = ISDestroy(&ISlocal); CHKERRV(ierr);
  ierr = ISDestroy(&ISlocalInv); CHKERRV(ierr);
  ierr = ISDestroy(&ISglobal); CHKERRV(ierr);
  ierr = ISDestroy(&ISidentity); CHKERRV(ierr);
  ierr = ISDestroy(&ISidentityGlobal); CHKERRV(ierr);
  ierr = ISLocalToGlobalMappingDestroy(&ISidentityMapping); CHKERRV(ierr);
  ierr = ISDestroy(&ISpolyparams); CHKERRV(ierr);
  ierr = ISLocalToGlobalMappingDestroy(&ISpolyparamsMapping); CHKERRV(ierr);

  const PetscInt *mapIndizes;
  ISLocalToGlobalMappingGetIndices(_ISmapping, &mapIndizes);

  Eigen::VectorXd distance(dimensions);

  // We do preallocating of the matrices C and A. That means we traverse the input data once, just
  // to know where we have entries in the sparse matrix. This information petsc can use to
  // preallocate the matrix. In the second phase we actually fill the matrix.

  // Stores col -> value for each row;
  std::vector<std::vector<std::pair<int, double>>> vertexData;

  if (_preallocation == Preallocation::SAVE) {
    vertexData = savedPreallocationMatrixC(inMesh);
  }
  if (_preallocation == Preallocation::COMPUTE) {
    computePreallocationMatrixC(inMesh);
  }
  if (_preallocation == Preallocation::ESTIMATE) {
    estimatePreallocationMatrixC(n, n, inMesh);
  }
  if (_preallocation == Preallocation::TREE) {
    vertexData = bgPreallocationMatrixC(inMesh);
  }

  // -- BEGIN FILL LOOP FOR MATRIX C --
  precice::utils::Event eFillC("map.pet.fillC.From" + input()->getName() + "To"+ output()->getName(), precice::syncMode);
  // We collect entries for each row and set them blockwise using MatSetValues.
  int preallocRow = 0;
  for (const mesh::Vertex& inVertex : inMesh->vertices()) {
    if (not inVertex.isOwner())
      continue;

    PetscInt const row = inVertex.getGlobalIndex() + polyparams;
    PetscInt colNum = 0;  // holds the number of columns
    PetscInt colIdx[_matrixC.getSize().second];     // holds the columns indices of the entries
    PetscScalar rowVals[_matrixC.getSize().second]; // holds the values of the entries

    // -- SETS THE POLYNOMIAL PART OF THE MATRIX --
    if (_polynomial == Polynomial::ON or _polynomial == Polynomial::SEPARATE) {
      colIdx[colNum] = colNum;
      rowVals[colNum++] = 1;

      for (int dim = 0; dim < dimensions; dim++) {
        if (not _deadAxis[dim]) {
          colIdx[colNum] = colNum;
          rowVals[colNum++] = inVertex.getCoords()[dim];
        }
      }

      if (_polynomial == Polynomial::ON) {
        ierr = MatSetValuesLocal(_matrixC, colNum, colIdx, 1, &row, rowVals, INSERT_VALUES); CHKERRV(ierr);
      }
      else if (_polynomial == Polynomial::SEPARATE) {
        ierr = MatSetValuesLocal(_matrixQ, 1, &row, colNum, colIdx, rowVals, INSERT_VALUES); CHKERRV(ierr);
      }
      colNum = 0;
    }

    // -- SETS THE COEFFICIENTS --
    if (_preallocation == Preallocation::SAVE or _preallocation == Preallocation::TREE) {
      const auto & rowVertices = vertexData[preallocRow];
      for (const auto & vertex : rowVertices) {
        rowVals[colNum] = _basisFunction.evaluate(vertex.second);
        colIdx[colNum++] = vertex.first;
      }
      ++preallocRow;
    }
    else {
      for (const mesh::Vertex& vj : inMesh->vertices()) {
        int col = vj.getGlobalIndex() + polyparams;
        if (row > col)
          continue; // matrix is symmetric
        distance = inVertex.getCoords() - vj.getCoords();
        for (int d = 0; d < dimensions; d++) {
          if (_deadAxis[d]) {
            distance[d] = 0;
          }
        }
        if (_basisFunction.getSupportRadius() > distance.norm()) {
          double coeff = _basisFunction.evaluate(distance.norm());
          rowVals[colNum] = coeff;
          colIdx[colNum++] = col; // column of entry is the globalIndex
        }
      }
    }
    ierr = MatSetValuesLocal(_matrixC, 1, &row, colNum, colIdx, rowVals, INSERT_VALUES); CHKERRV(ierr);
  }
  DEBUG("Finished filling Matrix C");
  eFillC.stop();
  // -- END FILL LOOP FOR MATRIX C --

  // PETSc requires that all diagonal entries are set, even if set to zero.
  _matrixC.assemble(MAT_FLUSH_ASSEMBLY);
  petsc::Vector zeros(_matrixC);
  MatDiagonalSet(_matrixC, zeros, ADD_VALUES);

  // Begin assembly here, all assembly is ended at the end of this function.
  ierr = MatAssemblyBegin(_matrixC, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
  ierr = MatAssemblyBegin(_matrixQ, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);

  if (_preallocation == Preallocation::SAVE) {
    vertexData = savedPreallocationMatrixA(inMesh, outMesh);
  }
  if (_preallocation == Preallocation::COMPUTE) {
    computePreallocationMatrixA(inMesh, outMesh);
  }
  if (_preallocation == Preallocation::ESTIMATE) {
    estimatePreallocationMatrixA(outputSize, n, inMesh);
  }
  if (_preallocation == Preallocation::TREE) {
    vertexData = bgPreallocationMatrixA(inMesh, outMesh);
  }

  // -- BEGIN FILL LOOP FOR MATRIX A --
  DEBUG("Begin filling matrix A.");
  precice::utils::Event eFillA("map.pet.fillA.From" + input()->getName() + "To" + output()->getName(), precice::syncMode);

  for (PetscInt row = ownerRangeABegin; row < ownerRangeAEnd; row++) {
    PetscInt colNum = 0;
    PetscInt colIdx[_matrixA.getSize().second];     // holds the columns indices of the entries
    PetscScalar rowVals[_matrixA.getSize().second]; // holds the values of the entries
    const mesh::Vertex& oVertex = outMesh->vertices()[row - _matrixA.ownerRange().first];

    // -- SET THE POLYNOMIAL PART OF THE MATRIX --
    if (_polynomial == Polynomial::ON or _polynomial == Polynomial::SEPARATE) {
      Mat m = _polynomial == Polynomial::ON ? _matrixA : _matrixV;

      colIdx[colNum] = colNum;
      rowVals[colNum++] = 1;

      for (int dim = 0; dim < dimensions; dim++) {
        if (not _deadAxis[dim]) {
          colIdx[colNum] = colNum;
          rowVals[colNum++] = oVertex.getCoords()[dim];
        }
      }
      ierr = MatSetValuesLocal(m, 1, &row, colNum, colIdx, rowVals, INSERT_VALUES); CHKERRV(ierr);
      colNum = 0;
    }

    // -- SETS THE COEFFICIENTS --
    if (_preallocation == Preallocation::SAVE or _preallocation == Preallocation::TREE) {
      const auto & rowVertices = vertexData[row - ownerRangeABegin];
      for (const auto & vertex : rowVertices) {
        rowVals[colNum] = _basisFunction.evaluate(vertex.second);
        colIdx[colNum++] = vertex.first;
      }
    }
    else {
      for (const mesh::Vertex& inVertex : inMesh->vertices()) {
        distance = oVertex.getCoords() - inVertex.getCoords();
        for (int d = 0; d < dimensions; d++) {
          if (_deadAxis[d])
            distance[d] = 0;
        }
        if (_basisFunction.getSupportRadius() > distance.norm()) {
          double coeff = _basisFunction.evaluate(distance.norm());
          rowVals[colNum] = coeff;
          colIdx[colNum++] = inVertex.getGlobalIndex() + polyparams;
        }
      }
    }
    ierr = MatSetValuesLocal(_matrixA, 1, &row, colNum, colIdx, rowVals, INSERT_VALUES); CHKERRV(ierr);
  }
  DEBUG("Finished filling Matrix A");
  eFillA.stop();
  // -- END FILL LOOP FOR MATRIX A --

  precice::utils::Event ePostFill("map.pet.postFill.From" + input()->getName() + "To" + output()->getName(), precice::syncMode);

  ierr = MatAssemblyBegin(_matrixA, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);

  ISLocalToGlobalMappingRestoreIndices(_ISmapping, &mapIndizes);
  ierr = MatAssemblyEnd(_matrixC, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
  ierr = MatAssemblyEnd(_matrixQ, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
  ierr = MatAssemblyEnd(_matrixA, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
  _matrixQ.assemble();
  _matrixV.assemble();

  // -- CONFIGURE SOLVER FOR POLYNOMIAL --
  if (_polynomial == Polynomial::SEPARATE) {
    PC pc;
    KSPGetPC(_QRsolver, &pc);
    PCSetType(pc, PCNONE);
    KSPSetType(_QRsolver, KSPLSQR);
    KSPSetOperators(_QRsolver, _matrixQ, _matrixQ);
  }

  // -- CONFIGURE SOLVER FOR SYSTEM MATRIX --
  KSPSetOperators(_solver, _matrixC, _matrixC); CHKERRV(ierr);
  KSPSetTolerances(_solver, _solverRtol, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
  KSPSetInitialGuessNonzero(_solver, PETSC_TRUE); CHKERRV(ierr); // Reuse the results from the last iteration, held in the out vector.
  KSPSetFromOptions(_solver);

  // if (totalNNZ > static_cast<size_t>(20*n)) {
  //   DEBUG("Using Cholesky decomposition as direct solver for dense matrix.");
  //   PC prec;
  //   KSPSetType(_solver, KSPPREONLY);
  //   KSPGetPC(_solver, &prec);
  //   PCSetType(prec, PCCHOLESKY);
  //   PCFactorSetShiftType(prec, MAT_SHIFT_NONZERO);
  // }

  // -- COMPUTE RESCALING COEFFICIENTS USING THE SYSTEM MATRIX SOLVER --
  if (useRescaling and (_polynomial == Polynomial::SEPARATE)) {
    petsc::Vector rhs(_matrixC);
    ierr = MatCreateVecs(_matrixC, nullptr, &rescalingCoeffs.vector); CHKERRV(ierr);
    VecSet(rhs, 1);
    rhs.assemble();
    _solver.solve(rhs, rescalingCoeffs);
  }

  _hasComputedMapping = true;

  DEBUG("Number of mallocs for matrix C = " << _matrixC.getInfo(MAT_LOCAL).mallocs);
  DEBUG("Non-zeros allocated / used / unused for matrix C = " << _matrixC.getInfo(MAT_LOCAL).nz_allocated << " / " << _matrixC.getInfo(MAT_LOCAL).nz_used << " / " << _matrixC.getInfo(MAT_LOCAL).nz_unneeded);
  DEBUG("Number of mallocs for matrix A = " << _matrixA.getInfo(MAT_LOCAL).mallocs);
  DEBUG("Non-zeros allocated / used / unused for matrix A = " << _matrixA.getInfo(MAT_LOCAL).nz_allocated << " / " << _matrixA.getInfo(MAT_LOCAL).nz_used << " / " << _matrixA.getInfo(MAT_LOCAL).nz_unneeded);
}

template<typename RADIAL_BASIS_FUNCTION_T>
bool PetRadialBasisFctMapping<RADIAL_BASIS_FUNCTION_T>:: hasComputedMapping() const
{
  return _hasComputedMapping;
}

template<typename RADIAL_BASIS_FUNCTION_T>
void PetRadialBasisFctMapping<RADIAL_BASIS_FUNCTION_T>:: clear()
{
  _matrixC.reset();
  _matrixA.reset();
  _matrixQ.reset();
  _matrixV.reset();

  _solver.reset();
  _QRsolver.reset();

  petsc::destroy(&_ISmapping);

  previousSolution.clear();
  _hasComputedMapping = false;
}

template<typename RADIAL_BASIS_FUNCTION_T>
void PetRadialBasisFctMapping<RADIAL_BASIS_FUNCTION_T>::map(int inputDataID, int outputDataID)
{
  TRACE(inputDataID, outputDataID);
  precice::utils::Event e("map.pet.mapData.From" + input()->getName() + "To" + output()->getName(), precice::syncMode);

  assertion(_hasComputedMapping);
  assertion(input()->getDimensions() == output()->getDimensions(),
            input()->getDimensions(), output()->getDimensions());

  PetscErrorCode ierr = 0;
  auto& inValues = input()->data(inputDataID)->values();
  auto& outValues = output()->data(outputDataID)->values();

  int valueDim = input()->data(inputDataID)->getDimensions();
  assertion(valueDim == output()->data(outputDataID)->getDimensions(),
            valueDim, output()->data(outputDataID)->getDimensions());

  if (getConstraint() == CONSERVATIVE) {
    petsc::Vector au(_matrixA, "au", petsc::Vector::RIGHT);
    petsc::Vector in(_matrixA, "in");
    int inRangeStart, inRangeEnd;
    std::tie(inRangeStart, inRangeEnd) = in.ownerRange();
    for (int dim = 0; dim < valueDim; dim++) {
      printMappingInfo(inputDataID, dim);

      // Fill input from input data values
      for (size_t i = 0; i < input()->vertices().size(); i++ ) {
        int globalIndex = input()->vertices()[i].getGlobalIndex(); // globalIndex is target row
        if (globalIndex >= inRangeStart and globalIndex < inRangeEnd) // only fill local rows
          ierr = VecSetValue(in, globalIndex, inValues[i*valueDim + dim], INSERT_VALUES); CHKERRV(ierr);
      }
      in.assemble();

      // Gets the petsc::vector for the given combination of outputData, inputData and dimension
      // If none created yet, create one, based on _matrixC
      petsc::Vector& out = std::get<0>(
        previousSolution.emplace(std::piecewise_construct,
                                 std::forward_as_tuple(inputDataID + outputDataID * 10 + dim * 100),
                                 std::forward_as_tuple(_matrixC, "out"))
        )->second;

      if (_polynomial == Polynomial::SEPARATE) {
        petsc::Vector epsilon(_matrixV, "epsilon", petsc::Vector::RIGHT);
        ierr = MatMultTranspose(_matrixV, in, epsilon); CHKERRV(ierr);
        petsc::Vector eta(_matrixA, "eta", petsc::Vector::RIGHT);
        ierr = MatMultTranspose(_matrixA, in, eta); CHKERRV(ierr);
        petsc::Vector mu(_matrixC, "mu", petsc::Vector::LEFT);
        _solver.solve(eta, mu);
        VecScale(epsilon, -1);
        petsc::Vector tau(_matrixQ, "tau", petsc::Vector::RIGHT);
        ierr = MatMultTransposeAdd(_matrixQ, mu, epsilon, tau); CHKERRV(ierr);
        petsc::Vector sigma(_matrixQ, "sigma", petsc::Vector::LEFT);
        ierr = KSPSolveTranspose(_QRsolver, tau, sigma); CHKERRV(ierr);
        VecWAXPY(out, -1, sigma, mu);
      }
      else {
        ierr = MatMultTranspose(_matrixA, in, au); CHKERRV(ierr);
        utils::Event eSolve("map.pet.solveConservative.From" + input()->getName() + "To" + output()->getName(), precice::syncMode);
        if (not _solver.solve(au, out)) {
          KSPView(_solver, PETSC_VIEWER_STDOUT_WORLD);
          ERROR("RBF linear system has not converged.");
        }
        eSolve.stop();

      }
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
    petsc::Vector a(_matrixQ, "a", petsc::Vector::RIGHT); // holds the solution of the LS polynomial

    ierr = VecSetLocalToGlobalMapping(in, _ISmapping); CHKERRV(ierr);
    const PetscScalar *vecArray;

    // For every data dimension, perform mapping
    for (int dim=0; dim < valueDim; dim++) {
      printMappingInfo(inputDataID, dim);

      // Fill input from input data values
      int count = 0;
      for (const auto& vertex : input()->vertices()) {
        ierr = VecSetValueLocal(in, vertex.getGlobalIndex()+polyparams, inValues[count*valueDim + dim], INSERT_VALUES); CHKERRV(ierr);
        count++;
      }
      in.assemble();

      if (_polynomial == Polynomial::SEPARATE) {
        if (not _QRsolver.solve(in, a)) {
          KSPView(_QRsolver, PETSC_VIEWER_STDOUT_WORLD);
          ERROR("Polynomial QR linear system has not converged.");
        }
        VecScale(a, -1);
        MatMultAdd(_matrixQ, a, in, in); // Subtract the polynomial from the input values
      }

      petsc::Vector& p = std::get<0>(  // Save and reuse the solution from the previous iteration
        previousSolution.emplace(std::piecewise_construct,
                                 std::forward_as_tuple(inputDataID + outputDataID * 10 + dim * 100),
                                 std::forward_as_tuple(_matrixC, "p"))
        )->second;

      utils::Event eSolve("map.pet.solveConsistent.From" + input()->getName() + "To" + output()->getName(), precice::syncMode);
      if (not _solver.solve(in, p)) {
        KSPView(_solver, PETSC_VIEWER_STDOUT_WORLD);
        ERROR("RBF linear system has not converged.");
      }
      eSolve.stop();

      ierr = MatMult(_matrixA, p, out); CHKERRV(ierr);

      if (useRescaling and (_polynomial == Polynomial::SEPARATE)) {
        petsc::Vector temp(_matrixA);
        PetscReal min;
        ierr = MatMult(_matrixA, rescalingCoeffs, temp); CHKERRV(ierr); // get the output of g(x) = 1
        ierr = VecMin(temp, nullptr, &min); CHKERRV(ierr); // check for zeros before devision
        if (min < 1e-6)
          WARN("Zero elements found in rescaling, omit rescaling.");
        else
          ierr = VecPointwiseDivide(out, out, temp); CHKERRV(ierr);
      }

      if (_polynomial == Polynomial::SEPARATE) {
        ierr = VecScale(a, -1); // scale it back to add the polynomial
        ierr = MatMultAdd(_matrixV, a, out, out); CHKERRV(ierr);
      }
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

/*
 * For the re-partitioning process with RBF mappings, also compare Figure 69 in Benjamin U's thesis (page 89).
 * https://mediatum.ub.tum.de/doc/1320661/document.pdf
 */
template<typename RADIAL_BASIS_FUNCTION_T>
void PetRadialBasisFctMapping<RADIAL_BASIS_FUNCTION_T>::tagMeshFirstRound()
{
  TRACE();
  mesh::PtrMesh filterMesh, otherMesh;
  if (getConstraint() == CONSISTENT){
    filterMesh = input(); // remote
    otherMesh = output(); // local
  }
  else if (getConstraint() == CONSERVATIVE) {
    filterMesh = output(); // remote
    otherMesh = input(); // local
  }

  if (otherMesh->vertices().size() == 0)
      return; // Ranks not at the interface should never hold interface vertices

  // Tags all vertices that are inside otherMesh's bounding box, enlarged by the support radius
  for(mesh::Vertex& v : filterMesh->vertices()) {
    bool isInside = true;
    #if PETSC_MAJOR >= 3 and PETSC_MINOR >= 8
    if (_basisFunction.hasCompactSupport()) {
      for (int d = 0; d < v.getDimensions(); d++) {
        if (v.getCoords()[d] < otherMesh->getBoundingBox()[d].first - _basisFunction.getSupportRadius() or
            v.getCoords()[d] > otherMesh->getBoundingBox()[d].second + _basisFunction.getSupportRadius() ) {
          isInside = false;
          break;
        }
      }
    }
    #else
      #warning "Mesh filtering deactivated, due to PETSc version < 3.8. \
preCICE is fully functional, but performance for large cases is degraded."
    #endif
    if (isInside)
      v.tag();
  }
}


/*
 * For the re-partitioning process with RBF mappings, also compare Figure 69 in Benjamin U's thesis (page 89).
 * https://mediatum.ub.tum.de/doc/1320661/document.pdf
 */
template<typename RADIAL_BASIS_FUNCTION_T>
void PetRadialBasisFctMapping<RADIAL_BASIS_FUNCTION_T>::tagMeshSecondRound()
{
  TRACE();

  if (not _basisFunction.hasCompactSupport())
    return; // Tags should not be changed

  mesh::PtrMesh mesh; // The mesh we want to filter

  if (getConstraint() == CONSISTENT)
    mesh = input();
  else if (getConstraint() == CONSERVATIVE)
    mesh = output();

  mesh::Mesh::BoundingBox bb(mesh->getDimensions(),
                             std::make_pair(std::numeric_limits<double>::max(),
                                            std::numeric_limits<double>::lowest()));

  // Construct bounding box around all owned vertices
  for (mesh::Vertex& v : mesh->vertices()) {
    if (v.isOwner()) {
      assertion(v.isTagged()); // Should be tagged from the first round
      for (int d = 0; d < v.getDimensions(); d++) {
        bb[d].first  = std::min(v.getCoords()[d], bb[d].first);
        bb[d].second = std::max(v.getCoords()[d], bb[d].second);
      }
    }
  }
  // Tag vertices that are inside the bounding box + support radius
  for (mesh::Vertex& v : mesh->vertices()) {
    bool isInside = true;
    for (int d = 0; d < v.getDimensions(); d++) {
      if (v.getCoords()[d] < bb[d].first - _basisFunction.getSupportRadius() or
          v.getCoords()[d] > bb[d].second + _basisFunction.getSupportRadius() ) {
        isInside = false;
        break;
      }
    }
    if (isInside)
      v.tag();
  }
}


template <typename RADIAL_BASIS_FUNCTION_T>
void PetRadialBasisFctMapping<RADIAL_BASIS_FUNCTION_T>::printMappingInfo(int inputDataID, int dim) const
{
  if (not precice::testMode) {
    const std::string constraintName = getConstraint() == CONSERVATIVE ? "conservative" : "consistent";
    const std::string polynomialName = _polynomial == Polynomial::ON ? "on" : _polynomial == Polynomial::OFF ? "off" : "separate";

    INFO("Mapping " << input()->data(inputDataID)->getName() << " " << constraintName
         << " from " << input()->getName() << " (ID " << input()->getID() << ")"
         << " to " << output()->getName() << " (ID " << output()->getID() << ") "
         << "for dimension " << dim << ") with polynomial set to " << polynomialName);
  }
}

template <typename RADIAL_BASIS_FUNCTION_T>
void PetRadialBasisFctMapping<RADIAL_BASIS_FUNCTION_T>::estimatePreallocationMatrixC(
  int rows, int cols, mesh::PtrMesh mesh)
{
  std::ignore = rows;
  std::ignore = cols;

  // auto rank = utils::Parallel::getProcessRank();
  auto size = utils::Parallel::getCommunicatorSize();
  auto supportRadius = _basisFunction.getSupportRadius();

  auto bbox = mesh->getBoundingBox();
  auto meshSize = mesh->vertices().size();

  double meshArea = 1;
  // WARN(bbox);
  for (int d = 0; d < getDimensions(); d++)
    if (not _deadAxis[d])
      meshArea *= bbox[d].second - bbox[d].first;

  // supportVolume = math::PI * 4.0/3.0 * std::pow(supportRadius, 3);
  double supportVolume = 0;
  if (getDimensions() == 2)
    supportVolume = 2 * supportRadius;
  else if (getDimensions() == 3)
    supportVolume = math::PI * std::pow(supportRadius, 2);

  WARN("supportVolume = " << supportVolume);
  WARN("meshArea = " << meshArea);
  WARN("meshSize = " << meshSize);

  int nnzPerRow = meshSize / meshArea * supportVolume;
  WARN("nnzPerRow = " << nnzPerRow);
  // int nnz = nnzPerRow * rows;

  if (utils::Parallel::getCommunicatorSize() == 1) {
    MatSeqSBAIJSetPreallocation(_matrixC, _matrixC.blockSize(), nnzPerRow / 2, nullptr);
  }
  else {
    MatMPISBAIJSetPreallocation(_matrixC, _matrixC.blockSize(),
                                nnzPerRow / (size * 2), nullptr , nnzPerRow / (size * 2), nullptr);
  }
  MatSetOption(_matrixC, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
}

template <typename RADIAL_BASIS_FUNCTION_T>
void PetRadialBasisFctMapping<RADIAL_BASIS_FUNCTION_T>::estimatePreallocationMatrixA(
  int rows, int cols, mesh::PtrMesh mesh)
{
  std::ignore = rows;
  std::ignore = cols;

  // auto rank = utils::Parallel::getProcessRank();
  auto size = utils::Parallel::getCommunicatorSize();
  auto supportRadius = _basisFunction.getSupportRadius();

  auto bbox = mesh->getBoundingBox();
  auto meshSize = mesh->vertices().size();

  double meshArea = 1;
  // WARN(bbox);
  for (int d = 0; d < getDimensions(); d++)
    if (not _deadAxis[d])
      meshArea *= bbox[d].second - bbox[d].first;

  // supportVolume = math::PI * 4.0/3.0 * std::pow(supportRadius, 3);
  double supportVolume = 0;
  if (getDimensions() == 2)
    supportVolume = 2 * supportRadius;
  else if (getDimensions() == 3)
    supportVolume = math::PI * std::pow(supportRadius, 2);

  WARN("supportVolume = " << supportVolume);
  WARN("meshArea = " << meshArea);
  WARN("meshSize = " << meshSize);

  int nnzPerRow = meshSize / meshArea * supportVolume;
  WARN("nnzPerRow = " << nnzPerRow);
  // int nnz = nnzPerRow * rows;

  if (utils::Parallel::getCommunicatorSize() == 1) {
    MatSeqSBAIJSetPreallocation(_matrixA, _matrixA.blockSize(), nnzPerRow, nullptr);
  }
  else {
    MatMPISBAIJSetPreallocation(_matrixA, _matrixA.blockSize(),
                                nnzPerRow / size, nullptr , nnzPerRow / size, nullptr);
  }
  MatSetOption(_matrixA, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
}


template <typename RADIAL_BASIS_FUNCTION_T>
void PetRadialBasisFctMapping<RADIAL_BASIS_FUNCTION_T>::computePreallocationMatrixC(const mesh::PtrMesh inMesh)
{
  precice::utils::Event ePreallocC("map.pet.preallocC.From" + input()->getName() + "To" + output()->getName(), precice::syncMode);

  PetscInt n, ierr;

  const PetscInt *mapIndizes;
  ISLocalToGlobalMappingGetIndices(_ISmapping, &mapIndizes);

  int dimensions = input()->getDimensions();
  Eigen::VectorXd distance(dimensions);

  std::tie(n, std::ignore) = _matrixC.getLocalSize();
  std::vector<PetscInt> d_nnz(n), o_nnz(n);
  PetscInt colOwnerRangeCBegin, colOwnerRangeCEnd;
  std::tie(colOwnerRangeCBegin, colOwnerRangeCEnd) = _matrixC.ownerRangeColumn();

  size_t local_row = 0;
  // -- PREALLOCATES THE POLYNOMIAL PART OF THE MATRIX --
  if (_polynomial == Polynomial::ON) {
    for (local_row = 0; local_row < localPolyparams; local_row++) {
      d_nnz[local_row] = colOwnerRangeCEnd - colOwnerRangeCBegin;
      o_nnz[local_row] = _matrixC.getSize().first - d_nnz[local_row];
    }
  }

  for (const mesh::Vertex& inVertex : inMesh->vertices()) {
    if (not inVertex.isOwner())
      continue;

    PetscInt col = polyparams - 1;
    const int global_row = local_row + _matrixC.ownerRange().first;
    d_nnz[local_row] = 0;
    o_nnz[local_row] = 0;

    // -- PREALLOCATES THE COEFFICIENTS --
    for (mesh::Vertex& vj : inMesh->vertices()) {
      col++;

      const int mapped_col = mapIndizes[vj.getGlobalIndex() + polyparams];
      if (global_row > mapped_col) // Skip, since we are below the diagonal
        continue;

      distance = inVertex.getCoords() - vj.getCoords();
      for (int d = 0; d < dimensions; d++) {
        if (_deadAxis[d]) {
          distance[d] = 0;
        }
      }

      if (_basisFunction.getSupportRadius() > distance.norm() or col == global_row) {
        if (mapped_col >= colOwnerRangeCBegin and mapped_col < colOwnerRangeCEnd)
          d_nnz[local_row]++;
        else
          o_nnz[local_row]++;
      }
    }
    local_row++;
  }

  if (utils::Parallel::getCommunicatorSize() == 1) {
    // std::cout << "Computed Preallocation C Seq diagonal = " << std::accumulate(d_nnz.begin(), d_nnz.end(), 0) << std::endl;
    ierr = MatSeqSBAIJSetPreallocation(_matrixC, _matrixC.blockSize(), 0, d_nnz.data()); CHKERRV(ierr);
  }
  else {
    // std::cout << "Computed Preallocation C MPI diagonal = " << std::accumulate(d_nnz.begin(), d_nnz.end(), 0)
    // << ", off-diagonal = " << std::accumulate(o_nnz.begin(), o_nnz.end(), 0) << std::endl;
    ierr = MatMPISBAIJSetPreallocation(_matrixC, _matrixC.blockSize(), 0, d_nnz.data(), 0, o_nnz.data()); CHKERRV(ierr);
  }
  MatSetOption(_matrixC, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);

  ISLocalToGlobalMappingRestoreIndices(_ISmapping, &mapIndizes);

  ePreallocC.stop();
}

template <typename RADIAL_BASIS_FUNCTION_T>
void PetRadialBasisFctMapping<RADIAL_BASIS_FUNCTION_T>::computePreallocationMatrixA(
  const mesh::PtrMesh inMesh, const mesh::PtrMesh outMesh)
{
  precice::utils::Event ePreallocA("map.pet.preallocA.From" + input()->getName() + "To" + output()->getName(), precice::syncMode);

  PetscInt ownerRangeABegin, ownerRangeAEnd, colOwnerRangeABegin, colOwnerRangeAEnd;
  PetscInt outputSize, ierr;

  const PetscInt *mapIndizes;
  ISLocalToGlobalMappingGetIndices(_ISmapping, &mapIndizes);

  std::tie(outputSize, std::ignore) = _matrixA.getLocalSize();

  std::vector<PetscInt> d_nnz(outputSize), o_nnz(outputSize);

  std::tie(ownerRangeABegin, ownerRangeAEnd) = _matrixA.ownerRange();
  std::tie(colOwnerRangeABegin, colOwnerRangeAEnd) = _matrixA.ownerRangeColumn();

  int dimensions = input()->getDimensions();
  Eigen::VectorXd distance(dimensions);

  for (int localRow = 0; localRow < ownerRangeAEnd - ownerRangeABegin; localRow++) {
    d_nnz[localRow] = 0;
    o_nnz[localRow] = 0;
    PetscInt col = 0;
    const mesh::Vertex& oVertex = outMesh->vertices()[localRow];

    // -- PREALLOCATE THE POLYNOM PART OF THE MATRIX --
    if (_polynomial == Polynomial::ON) {
      if (mapIndizes[col] >= colOwnerRangeABegin and mapIndizes[col] < colOwnerRangeAEnd)
        d_nnz[localRow]++;
      else
        o_nnz[localRow]++;
      col++;

      for (int dim = 0; dim < dimensions; dim++) {
        if (not _deadAxis[dim]) {
          if (mapIndizes[col] >= colOwnerRangeABegin and mapIndizes[col] < colOwnerRangeAEnd)
            d_nnz[localRow]++;
          else
            o_nnz[localRow]++;
          col++;
        }
      }
    }

    // -- PREALLOCATE THE COEFFICIENTS --
    for (const mesh::Vertex& inVertex : inMesh->vertices()) {
      distance = oVertex.getCoords() - inVertex.getCoords();
      for (int d = 0; d < dimensions; d++) {
        if (_deadAxis[d])
          distance[d] = 0;
      }

      if (_basisFunction.getSupportRadius() > distance.norm()) {
        col = inVertex.getGlobalIndex() + polyparams;
        if (mapIndizes[col] >= colOwnerRangeABegin and mapIndizes[col] < colOwnerRangeAEnd)
          d_nnz[localRow]++;
        else
          o_nnz[localRow]++;
      }
    }
  }
  if (utils::Parallel::getCommunicatorSize() == 1) {
    // std::cout << "Preallocation A Seq diagonal = " << std::accumulate(d_nnz.begin(), d_nnz.end(), 0) << std::endl;
    ierr = MatSeqAIJSetPreallocation(_matrixA, 0, d_nnz.data()); CHKERRV(ierr);
  }
  else {
    // std::cout << "Preallocation A MPI diagonal = " << std::accumulate(d_nnz.begin(), d_nnz.end(), 0)
    // << ", off-diagonal = " << std::accumulate(o_nnz.begin(), o_nnz.end(), 0) << std::endl;
    ierr = MatMPIAIJSetPreallocation(_matrixA, 0, d_nnz.data(), 0, o_nnz.data()); CHKERRV(ierr);
  }
  MatSetOption(_matrixA, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);

  ISLocalToGlobalMappingRestoreIndices(_ISmapping, &mapIndizes);

  ePreallocA.stop();
}

template <typename RADIAL_BASIS_FUNCTION_T>
std::vector<std::vector<std::pair<int, double>>> PetRadialBasisFctMapping<RADIAL_BASIS_FUNCTION_T>::savedPreallocationMatrixC(
  const mesh::PtrMesh inMesh)
{
  precice::utils::Event ePreallocC("map.pet.preallocC.From" + input()->getName() + "To" + output()->getName(), precice::syncMode);

  PetscInt n;

  const PetscInt *mapIndizes;
  ISLocalToGlobalMappingGetIndices(_ISmapping, &mapIndizes);

  int dimensions = input()->getDimensions();
  Eigen::VectorXd distance(dimensions);

  std::tie(n, std::ignore) = _matrixC.getLocalSize();
  std::vector<PetscInt> d_nnz(n), o_nnz(n);
  PetscInt colOwnerRangeCBegin, colOwnerRangeCEnd;
  std::tie(colOwnerRangeCBegin, colOwnerRangeCEnd) = _matrixC.ownerRangeColumn();

  std::vector<std::vector<std::pair<int, double>>> vertexData(n - localPolyparams);

  size_t local_row = 0;
  // -- PREALLOCATES THE POLYNOMIAL PART OF THE MATRIX --
  if (_polynomial == Polynomial::ON) {
    for (local_row = 0; local_row < localPolyparams; local_row++) {
      d_nnz[local_row] = colOwnerRangeCEnd - colOwnerRangeCBegin;
      o_nnz[local_row] = _matrixC.getSize().first - d_nnz[local_row];
    }
  }

  for (const mesh::Vertex& inVertex : inMesh->vertices()) {
    if (not inVertex.isOwner())
      continue;

    PetscInt col = polyparams - 1;
    const int global_row = local_row + _matrixC.ownerRange().first;
    d_nnz[local_row] = 0;
    o_nnz[local_row] = 0;

    // -- PREALLOCATES THE COEFFICIENTS --
    for (mesh::Vertex& vj : inMesh->vertices()) {
      col++;

      const int mapped_col = mapIndizes[vj.getGlobalIndex() + polyparams];
      if (global_row > mapped_col) // Skip, since we are below the diagonal
        continue;

      distance = inVertex.getCoords() - vj.getCoords();
      for (int d = 0; d < dimensions; d++)
        if (_deadAxis[d])
          distance[d] = 0;

      if (_basisFunction.getSupportRadius() > distance.norm() or col == global_row) {
        vertexData[local_row - localPolyparams].emplace_back( vj.getGlobalIndex() + polyparams, distance.norm() );
        if (mapped_col >= colOwnerRangeCBegin and mapped_col < colOwnerRangeCEnd)
          d_nnz[local_row]++;
        else
          o_nnz[local_row]++;
      }
    }
    local_row++;
  }

  if (utils::Parallel::getCommunicatorSize() == 1) {
    // std::cout << "Computed Preallocation C Seq diagonal = " << std::accumulate(d_nnz.begin(), d_nnz.end(), 0) << std::endl;
    MatSeqSBAIJSetPreallocation(_matrixC, _matrixC.blockSize(), 0, d_nnz.data());
  }
  else {
    // std::cout << "Computed Preallocation C MPI diagonal = " << std::accumulate(d_nnz.begin(), d_nnz.end(), 0)
    //           << ", off-diagonal = " << std::accumulate(o_nnz.begin(), o_nnz.end(), 0) << std::endl;
    MatMPISBAIJSetPreallocation(_matrixC, _matrixC.blockSize(), 0, d_nnz.data(), 0, o_nnz.data());
  }
  MatSetOption(_matrixC, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);
  ISLocalToGlobalMappingRestoreIndices(_ISmapping, &mapIndizes);

  ePreallocC.stop();

  return vertexData;
}

template <typename RADIAL_BASIS_FUNCTION_T>
std::vector<std::vector<std::pair<int, double>>> PetRadialBasisFctMapping<RADIAL_BASIS_FUNCTION_T>::savedPreallocationMatrixA(
  const mesh::PtrMesh inMesh, const mesh::PtrMesh outMesh)
{
  INFO("Using saved preallocation");
  precice::utils::Event ePreallocA("map.pet.preallocA.From" + input()->getName() + "To" + output()->getName(), precice::syncMode);

  PetscInt ownerRangeABegin, ownerRangeAEnd, colOwnerRangeABegin, colOwnerRangeAEnd;
  PetscInt outputSize, n;

  const PetscInt *mapIndizes;
  ISLocalToGlobalMappingGetIndices(_ISmapping, &mapIndizes);

  std::tie(outputSize, n) = _matrixA.getLocalSize(); // n hier nehmen
  std::tie(ownerRangeABegin, ownerRangeAEnd) = _matrixA.ownerRange();
  std::tie(colOwnerRangeABegin, colOwnerRangeAEnd) = _matrixA.ownerRangeColumn();
  int dimensions = input()->getDimensions();

  std::vector<PetscInt> d_nnz(outputSize), o_nnz(outputSize);
  Eigen::VectorXd distance(dimensions);

  // Contains localRow<localCols<colPosition, distance>>>
  std::vector<std::vector<std::pair<int, double>>> vertexData(outputSize);

  for (int localRow = 0; localRow < ownerRangeAEnd - ownerRangeABegin; localRow++) {
    d_nnz[localRow] = 0;
    o_nnz[localRow] = 0;
    PetscInt col = 0;
    const mesh::Vertex& oVertex = outMesh->vertices()[localRow];

    // -- PREALLOCATE THE POLYNOM PART OF THE MATRIX --
    if (_polynomial == Polynomial::ON) {
      if (mapIndizes[col] >= colOwnerRangeABegin and mapIndizes[col] < colOwnerRangeAEnd)
        d_nnz[localRow]++;
      else
        o_nnz[localRow]++;
      col++;

      for (int dim = 0; dim < dimensions; dim++) {
        if (not _deadAxis[dim]) {
          if (mapIndizes[col] >= colOwnerRangeABegin and mapIndizes[col] < colOwnerRangeAEnd)
            d_nnz[localRow]++;
          else
            o_nnz[localRow]++;
          col++;
        }
      }
    }

    // -- PREALLOCATE THE COEFFICIENTS --
    for (const mesh::Vertex& inVertex : inMesh->vertices()) {
      distance = oVertex.getCoords() - inVertex.getCoords();
      for (int d = 0; d < dimensions; d++)
        if (_deadAxis[d])
          distance[d] = 0;

      if (_basisFunction.getSupportRadius() > distance.norm()) {
        col = inVertex.getGlobalIndex() + polyparams;
        vertexData[localRow].emplace_back(col, distance.norm());

        if (mapIndizes[col] >= colOwnerRangeABegin and mapIndizes[col] < colOwnerRangeAEnd)
          d_nnz[localRow]++;
        else
          o_nnz[localRow]++;
      }
    }
  }
  if (utils::Parallel::getCommunicatorSize() == 1) {
    // std::cout << "Preallocation A Seq diagonal = " << std::accumulate(d_nnz.begin(), d_nnz.end(), 0) << std::endl;
    MatSeqAIJSetPreallocation(_matrixA, 0, d_nnz.data());
  }
  else {
    // std::cout << "Preallocation A MPI diagonal = " << std::accumulate(d_nnz.begin(), d_nnz.end(), 0)
    //           << ", off-diagonal = " << std::accumulate(o_nnz.begin(), o_nnz.end(), 0) << std::endl;
    MatMPIAIJSetPreallocation(_matrixA, 0, d_nnz.data(), 0, o_nnz.data());
  }
  MatSetOption(_matrixA, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);

  ISLocalToGlobalMappingRestoreIndices(_ISmapping, &mapIndizes);
  ePreallocA.stop();
  return vertexData;
}


template <typename RADIAL_BASIS_FUNCTION_T>
std::vector<std::vector<std::pair<int, double>>> PetRadialBasisFctMapping<RADIAL_BASIS_FUNCTION_T>::bgPreallocationMatrixC(
  const mesh::PtrMesh inMesh)
{
  INFO("Using tree-based preallocation for matrix C");
  precice::utils::Event ePreallocC("map.pet.preallocC.From" + input()->getName() + "To" + output()->getName(), precice::syncMode);
  namespace bg = boost::geometry;

  PetscInt n;

  auto tree = mesh::rtree::getVertexRTree(inMesh);
  double supportRadius = _basisFunction.getSupportRadius();

  const PetscInt *mapIndizes;
  ISLocalToGlobalMappingGetIndices(_ISmapping, &mapIndizes);

  int dimensions = input()->getDimensions();
  Eigen::VectorXd distance(dimensions);

  std::tie(n, std::ignore) = _matrixC.getLocalSize();
  std::vector<PetscInt> d_nnz(n), o_nnz(n);
  PetscInt colOwnerRangeCBegin, colOwnerRangeCEnd;
  std::tie(colOwnerRangeCBegin, colOwnerRangeCEnd) = _matrixC.ownerRangeColumn();

  std::vector<std::vector<std::pair<int, double>>> vertexData(n - localPolyparams);

  size_t local_row = 0;
  // -- PREALLOCATES THE POLYNOMIAL PART OF THE MATRIX --
  if (_polynomial == Polynomial::ON) {
    for (local_row = 0; local_row < localPolyparams; local_row++) {
      d_nnz[local_row] = colOwnerRangeCEnd - colOwnerRangeCBegin;
      o_nnz[local_row] = _matrixC.getSize().first - d_nnz[local_row];
    }
  }

  for (const mesh::Vertex& inVertex : inMesh->vertices()) {
    if (not inVertex.isOwner())
      continue;

    PetscInt col = polyparams - 1;
    const int global_row = local_row + _matrixC.ownerRange().first;
    d_nnz[local_row] = 0;
    o_nnz[local_row] = 0;

    // -- PREALLOCATES THE COEFFICIENTS --
    std::vector<size_t> results;
    auto search_box = mesh::getEnclosingBox(inVertex, supportRadius);

    tree->query(bg::index::within(search_box) and bg::index::satisfies([&](size_t const i){
          return bg::distance(inVertex, inMesh->vertices()[i]) <= supportRadius;}),
      std::back_inserter(results));

    // for (mesh::Vertex& vj : inMesh->vertices()) {
    for (auto i : results) {
      const mesh::Vertex & vj = inMesh->vertices()[i];
      col++;

      const int mapped_col = mapIndizes[vj.getGlobalIndex() + polyparams];
      if (global_row > mapped_col) // Skip, since we are below the diagonal
        continue;

      distance = inVertex.getCoords() - vj.getCoords();
      for (int d = 0; d < dimensions; d++)
        if (_deadAxis[d])
          distance[d] = 0;

      if (_basisFunction.getSupportRadius() > distance.norm() or col == global_row) {
        vertexData[local_row - localPolyparams].emplace_back( vj.getGlobalIndex() + polyparams, distance.norm() );
        if (mapped_col >= colOwnerRangeCBegin and mapped_col < colOwnerRangeCEnd)
          d_nnz[local_row]++;
        else
          o_nnz[local_row]++;
      }
    }
    local_row++;
  }

  if (utils::Parallel::getCommunicatorSize() == 1) {
    MatSeqSBAIJSetPreallocation(_matrixC, _matrixC.blockSize(), 0, d_nnz.data());
  }
  else {
    MatMPISBAIJSetPreallocation(_matrixC, _matrixC.blockSize(), 0, d_nnz.data(), 0, o_nnz.data());
  }
  MatSetOption(_matrixC, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);
  ISLocalToGlobalMappingRestoreIndices(_ISmapping, &mapIndizes);

  ePreallocC.stop();

  return vertexData;
}

template <typename RADIAL_BASIS_FUNCTION_T>
std::vector<std::vector<std::pair<int, double>>> PetRadialBasisFctMapping<RADIAL_BASIS_FUNCTION_T>::bgPreallocationMatrixA(
  const mesh::PtrMesh inMesh, const mesh::PtrMesh outMesh)
{
  INFO("Using tree-based preallocation for matrix A");
  precice::utils::Event ePreallocA("map.pet.preallocA.From" + input()->getName() + "To" + output()->getName(), precice::syncMode);
  namespace bg = boost::geometry;

  PetscInt ownerRangeABegin, ownerRangeAEnd, colOwnerRangeABegin, colOwnerRangeAEnd;
  PetscInt outputSize, n;
  auto tree = mesh::rtree::getVertexRTree(inMesh);
  double supportRadius = _basisFunction.getSupportRadius();

  const PetscInt *mapIndizes;
  ISLocalToGlobalMappingGetIndices(_ISmapping, &mapIndizes);

  std::tie(outputSize, n) = _matrixA.getLocalSize(); // n hier nehmen
  std::tie(ownerRangeABegin, ownerRangeAEnd) = _matrixA.ownerRange();
  std::tie(colOwnerRangeABegin, colOwnerRangeAEnd) = _matrixA.ownerRangeColumn();
  int dimensions = input()->getDimensions();

  std::vector<PetscInt> d_nnz(outputSize), o_nnz(outputSize);
  Eigen::VectorXd distance(dimensions);

  // Contains localRow<localCols<colPosition, distance>>>
  std::vector<std::vector<std::pair<int, double>>> vertexData(outputSize);

  for (int localRow = 0; localRow < ownerRangeAEnd - ownerRangeABegin; localRow++) {
    d_nnz[localRow] = 0;
    o_nnz[localRow] = 0;
    PetscInt col = 0;
    const mesh::Vertex& oVertex = outMesh->vertices()[localRow];

    // -- PREALLOCATE THE POLYNOM PART OF THE MATRIX --
    if (_polynomial == Polynomial::ON) {
      if (mapIndizes[col] >= colOwnerRangeABegin and mapIndizes[col] < colOwnerRangeAEnd)
        d_nnz[localRow]++;
      else
        o_nnz[localRow]++;
      col++;

      for (int dim = 0; dim < dimensions; dim++) {
        if (not _deadAxis[dim]) {
          if (mapIndizes[col] >= colOwnerRangeABegin and mapIndizes[col] < colOwnerRangeAEnd)
            d_nnz[localRow]++;
          else
            o_nnz[localRow]++;
          col++;
        }
      }
    }

    // -- PREALLOCATE THE COEFFICIENTS --
    std::vector<size_t> results;
    auto search_box = mesh::getEnclosingBox(oVertex, supportRadius);

    tree->query(bg::index::within(search_box) and bg::index::satisfies([&](size_t const i){
          return bg::distance(oVertex, inMesh->vertices()[i]) <= supportRadius;}),
      std::back_inserter(results));

    for (auto i : results ) {
        const mesh::Vertex & inVertex = inMesh->vertices()[i];
        distance = oVertex.getCoords() - inVertex.getCoords();

        for (int d = 0; d < dimensions; d++)
          if (_deadAxis[d])
            distance[d] = 0;

        if (_basisFunction.getSupportRadius() > distance.norm()) {
          col = inVertex.getGlobalIndex() + polyparams;
          vertexData[localRow].emplace_back(col, distance.norm());

          if (mapIndizes[col] >= colOwnerRangeABegin and mapIndizes[col] < colOwnerRangeAEnd)
            d_nnz[localRow]++;
          else
            o_nnz[localRow]++;
        }
    }
  }
  if (utils::Parallel::getCommunicatorSize() == 1) {
    MatSeqAIJSetPreallocation(_matrixA, 0, d_nnz.data());
  }
  else {
    MatMPIAIJSetPreallocation(_matrixA, 0, d_nnz.data(), 0, o_nnz.data());
  }
  MatSetOption(_matrixA, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);

  ISLocalToGlobalMappingRestoreIndices(_ISmapping, &mapIndizes);
  ePreallocA.stop();
  return vertexData;
}


}} // namespace precice, mapping

#endif // PRECICE_NO_PETSC
