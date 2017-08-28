#pragma once
#ifndef PRECICE_NO_PETSC

#include "mapping/Mapping.hpp"

#include <map>

#include "math/math.hpp"
#include "impl/BasisFunctions.hpp"
#include "config/MappingConfiguration.hpp"
#include "utils/Petsc.hpp"
namespace petsc = precice::utils::petsc;
#include "utils/EventTimings.hpp"

#include "petscmat.h"
#include "petscksp.h"
#include "petsclog.h"

// Forward declaration to friend the boost test struct
namespace MappingTests {
namespace PetRadialBasisFunctionMapping {
namespace Serial {
struct SolutionCaching;
}}}

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
   * @param[in] dimension Dimensionality of the meshes
   * @param[in] function Radial basis function used for mapping.
   * @param[in] xDead, yDead, zDead Deactivates mapping along an axis
   * @param[in] solverRtol Relative tolerance for the linear solver.
   * @param[in] polynomial Type of polynomial augmentation
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
    Polynomial              polynomial = Polynomial::ON);

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

  /// Logging device.
  static logging::Logger _log;

  bool _hasComputedMapping;

  /// Radial basis function type used in interpolation.
  RADIAL_BASIS_FUNCTION_T _basisFunction;

  /// Interpolation system matrix. Evaluated basis function on the input mesh
  petsc::Matrix _matrixC;

  /// Vandermonde Matrix, constructed from vertices of the input mesh
  petsc::Matrix _matrixQ;

  /// Interpolation evaluation matrix. Evaluated basis function on the output mesh
  petsc::Matrix _matrixA;

  /// Coordinates of the output mesh to evaluate the separated polynomial
  petsc::Matrix _matrixV;

  petsc::Vector rescalingCoeffs;

  KSP _solver;

  /// Used to solve the under-determined system for the separated polynomial.
  KSP _QRsolver;

  ISLocalToGlobalMapping _ISmapping;

  const double _solverRtol;

  /// true if the mapping along some axis should be ignored
  bool* _deadAxis;

  /// Toggles the use of the additonal polynomial
  const Polynomial _polynomial;

  /// Toggles use of rescaled basis functions, only active when Polynomial == SEPARATE
  const bool useRescaling = true;

  /// Number of coefficients for the integrated polynomial. Depends on dimension and number of dead dimensions
  size_t polyparams;

  /// Number of coefficients for the separated polynomial. Depends on dimension and number of dead dimensions
  size_t sepPolyparams;

  virtual bool doesVertexContribute(int vertexID) const override;

  void incPrealloc(PetscInt* diag, PetscInt* offDiag, int pos, int begin, int end);

  /// Caches the solution from the previous iteration, used as starting value for current iteration
  std::map<unsigned int, petsc::Vector> previousSolution;

  /// Prints an INFO about the current mapping
  void printMappingInfo(int inputDataID, int dim) const;
};

// --------------------------------------------------- HEADER IMPLEMENTATIONS

template<typename RADIAL_BASIS_FUNCTION_T>
logging::Logger PetRadialBasisFctMapping<RADIAL_BASIS_FUNCTION_T>::_log("mapping::PetRadialBasisFctMapping");

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
  Polynomial              polynomial)
  :
  Mapping ( constraint, dimensions ),
  _hasComputedMapping ( false ),
  _basisFunction ( function ),
  _matrixC(PETSC_COMM_WORLD, "C"),
  _matrixQ(PETSC_COMM_WORLD, "Q"),
  _matrixA(PETSC_COMM_WORLD, "A"),
  _matrixV(PETSC_COMM_WORLD, "V"),
  _solver(nullptr),
  _QRsolver(nullptr),
  _ISmapping(nullptr),
  _solverRtol(solverRtol),
  _polynomial(polynomial)
{
  setInputRequirement(VERTEX);
  setOutputRequirement(VERTEX);
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

  KSPCreate(PETSC_COMM_WORLD, &_solver);
  KSPCreate(PETSC_COMM_WORLD, &_QRsolver);
}

template<typename RADIAL_BASIS_FUNCTION_T>
PetRadialBasisFctMapping<RADIAL_BASIS_FUNCTION_T>::~PetRadialBasisFctMapping()
{
  delete[] _deadAxis;
  petsc::destroy(&_ISmapping);
  petsc::destroy(&_solver);
  petsc::destroy(&_QRsolver);
}


template<typename RADIAL_BASIS_FUNCTION_T>
void PetRadialBasisFctMapping<RADIAL_BASIS_FUNCTION_T>::computeMapping()
{
  TRACE();
  precice::utils::Event e(__func__);

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
  std::vector<int> myIndizes;

  // Indizes for Q^T, holding the polynom
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
  _matrixC.reset();
  _matrixC.init(n, n, PETSC_DETERMINE, PETSC_DETERMINE, MATSBAIJ);
  DEBUG("Set matrix C to local size " << n << " x " << n);
  ierr = MatSetOption(_matrixC, MAT_SYMMETRIC, PETSC_TRUE); CHKERRV(ierr);
  ierr = MatSetOption(_matrixC, MAT_SYMMETRY_ETERNAL, PETSC_TRUE); CHKERRV(ierr);

  // Matrix Q: Dense, holds the input mesh for the polynom if set to SEPERATE. Zero size otherwise
  _matrixQ.reset();
  _matrixQ.init(n, sepPolyparams, PETSC_DETERMINE, PETSC_DETERMINE, MATDENSE);
  DEBUG("Set matrix Q to local size " << n << " x " << sepPolyparams);

  // Matrix V: Dense, holds the output mesh for polynom if set to SEPERATE. Zero size otherwise
  _matrixV.reset();
  _matrixV.init(outputSize, sepPolyparams, PETSC_DETERMINE, PETSC_DETERMINE, MATDENSE);
  DEBUG("Set matrix V to local size " << outputSize << " x " << sepPolyparams);    
    
  // Matrix A: Sparse matrix with outputSize x n local size.
  _matrixA.reset();
  _matrixA.init(outputSize, n, PETSC_DETERMINE, PETSC_DETERMINE, MATAIJ);
  DEBUG("Set matrix A to local size " << outputSize << " x " << n);

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

  ierr = MatSetLocalToGlobalMapping(_matrixQ, _ISmapping, _ISmapping); CHKERRV(ierr);
  ierr = MatSetLocalToGlobalMapping(_matrixV, ISidentityMapping, _ISmapping); CHKERRV(ierr);
  
  // Destroy all local index sets and mappings
  ierr = ISDestroy(&ISlocal); CHKERRV(ierr);
  ierr = ISDestroy(&ISlocalInv); CHKERRV(ierr);
  ierr = ISDestroy(&ISglobal); CHKERRV(ierr);
  ierr = ISDestroy(&ISidentity); CHKERRV(ierr);
  ierr = ISDestroy(&ISidentityGlobal); CHKERRV(ierr);
  ierr = ISLocalToGlobalMappingDestroy(&ISidentityMapping); CHKERRV(ierr);

  Eigen::VectorXd distance(dimensions);

  // We do preallocating of the matrices C and A. That means we traverse the input data once, just
  // to know where we have entries in the sparse matrix. This information petsc can use to
  // preallocate the matrix. In the second phase we actually fill the matrix.

  // -- BEGIN PREALLOC LOOP FOR MATRIX C --
  // TODO Testen ob Preallocation perfekt ist.
  // MatSetOption(_matrixC, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  // MatSetOption(_matrixA, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  // {
  //   DEBUG("Begin preallocation matrix C");
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
  //     // DEBUG("Row = " << row << " n = " << n);
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
  //     // DEBUG("Reserved for row = " << row << ", reserved = " << nnz[row]);
  //   }
  //   PetscLogEventEnd(logPreallocCLoop, 0, 0, 0, 0);
  //   // -- END PREALLOC LOOP FOR MATRIX C --

  //   ierr = MatMPISBAIJSetPreallocation(_matrixC, 1, 0, nnz, 0, nnz); CHKERRV(ierr);    
  // }
  
  // -- BEGIN FILL LOOP FOR MATRIX C --
  int logCLoop = 2;
  PetscLogEventRegister("Filling Matrix C", 0, &logCLoop);
  PetscLogEventBegin(logCLoop, 0, 0, 0, 0);
  precice::utils::Event eFillC("PetRBF.fillC");
  // We collect entries for each row and set them blockwise using MatSetValues.
  for (const mesh::Vertex& inVertex : inMesh->vertices()) {
    if (not inVertex.isOwner())
      continue;

    int row = inVertex.getGlobalIndex() + polyparams;
    PetscInt colNum = 0;  // holds the number of columns
    PetscInt colIdx[_matrixC.getSize().second];     // holds the columns indices of the entries
    PetscScalar rowVals[_matrixC.getSize().second]; // holds the values of the entries
    
    // -- SETS THE POLYNOM PART OF THE MATRIX --
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
    for (mesh::Vertex& vj : inMesh->vertices()) {
      distance = inVertex.getCoords() - vj.getCoords();
      for (int d = 0; d < dimensions; d++) {
        if (_deadAxis[d]) {
          distance[d] = 0;
        }
      }
      double coeff = _basisFunction.evaluate(distance.norm());
      if (not math::equals(coeff, 0.0)) {
        rowVals[colNum] = coeff;
        colIdx[colNum] = vj.getGlobalIndex() + polyparams; // column of entry is the globalIndex
        colNum++;
      }
    }

    ierr = MatSetValuesLocal(_matrixC, 1, &row, colNum, colIdx, rowVals, INSERT_VALUES); CHKERRV(ierr);
  }
  DEBUG("Finished filling Matrix C");
  eFillC.stop();
  PetscLogEventEnd(logCLoop, 0, 0, 0, 0);
  // -- END FILL LOOP FOR MATRIX C --

  // Petsc requires that all diagonal entries are set, even if set to zero.
  _matrixC.assemble(MAT_FLUSH_ASSEMBLY);
  petsc::Vector zeros(_matrixC);
  MatDiagonalSet(_matrixC, zeros, ADD_VALUES);

  // Begin assembly here, all assembly is ended at the end of this function.
  ierr = MatAssemblyBegin(_matrixC, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
  ierr = MatAssemblyBegin(_matrixQ, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
  
  // -- BEGIN PREALLOC LOOP FOR MATRIX A --
  // {
  //   int localDiagColBegin = _matrixA.ownerRangeColumn().first;
  //   int localDiagColEnd = _matrixA.ownerRangeColumn().second;
  //   DEBUG("Local Submatrix Rows = " << ownerRangeABegin << " / " << ownerRangeAEnd <<
  //                ", Local Submatrix Cols = " << localDiagColBegin << " / " << localDiagColEnd);

  //   DEBUG("Begin preallocation matrix A.");
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
  //       DEBUG("Preallocating     diagonal row " << i << " with " << DnnzA[i] << " elements.");
  //       DEBUG("Preallocating off-diagonal row " << i << " with " << nnzA[i] << " elements.");
  //     }
  //   }
  //   // -- END PREALLOC LOOP FOR MATRIX A --

  //   // Preallocation: Number of diagonal non-zero entries equals outputSize, since diagonal is full
  //   MatSetOption(_matrixA, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  //   ierr = MatMPIAIJSetPreallocation(_matrixA, 0, DnnzA, 0, nnzA); CHKERRV(ierr);
  // }
  
  // -- BEGIN FILL LOOP FOR MATRIX A --
  DEBUG("Begin filling matrix A.");
  int logALoop = 4;
  PetscLogEventRegister("Filling Matrix A", 0, &logALoop);
  PetscLogEventBegin(logALoop, 0, 0, 0, 0);
  precice::utils::Event eFillA("PetRBF.fillA");

  for (int it = ownerRangeABegin; it < ownerRangeAEnd; it++) {
    PetscInt colNum = 0;
    PetscInt colIdx[_matrixA.getSize().second];     // holds the columns indices of the entries
    PetscScalar rowVals[_matrixA.getSize().second]; // holds the values of the entries
    const mesh::Vertex& oVertex = outMesh->vertices()[it - _matrixA.ownerRange().first];

    // -- SET THE POLYNOM PART OF THE MATRIX --
    if (_polynomial == Polynomial::ON or _polynomial == Polynomial::SEPARATE) {
      Mat m = _polynomial == Polynomial::ON ? _matrixA.matrix : _matrixV.matrix;
      
      colIdx[colNum] = colNum;
      rowVals[colNum++] = 1;
      
      for (int dim = 0; dim < dimensions; dim++) {
        if (not _deadAxis[dim]) {
          colIdx[colNum] = colNum;
          rowVals[colNum++] = oVertex.getCoords()[dim];
        }
      }
      ierr = MatSetValuesLocal(m, 1, &it, colNum, colIdx, rowVals, INSERT_VALUES); CHKERRV(ierr);
      colNum = 0;
    }
    
    // -- SETS THE COEFFICIENTS --
    for (const mesh::Vertex& inVertex : inMesh->vertices()) {
      distance = oVertex.getCoords() - inVertex.getCoords();
      for (int d = 0; d < dimensions; d++) {
        if (_deadAxis[d])
          distance[d] = 0;
      }
      double coeff = _basisFunction.evaluate(distance.norm());
      if (not math::equals(coeff, 0.0)) {
        rowVals[colNum] = coeff;
        colIdx[colNum] = inVertex.getGlobalIndex() + polyparams;
        colNum++;
      }
    }
    // DEBUG("Filling A: row = " << it << ", col count = " << colNum);
    ierr = MatSetValuesLocal(_matrixA, 1, &it, colNum, colIdx, rowVals, INSERT_VALUES); CHKERRV(ierr);
  }
  DEBUG("Finished filling Matrix A");
  eFillA.stop();
  PetscLogEventEnd(logALoop, 0, 0, 0, 0);
  // -- END FILL LOOP FOR MATRIX A --
  
  ierr = MatAssemblyBegin(_matrixA, MAT_FINAL_ASSEMBLY); CHKERRV(ierr);
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
    ierr = KSPSolve(_solver, rhs, rescalingCoeffs); CHKERRV(ierr);
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
  previousSolution.clear();
  _hasComputedMapping = false;
  // ISmapping destroy??
}

template<typename RADIAL_BASIS_FUNCTION_T>
void PetRadialBasisFctMapping<RADIAL_BASIS_FUNCTION_T>::map(int inputDataID, int outputDataID)
{
  TRACE(inputDataID, outputDataID);
  precice::utils::Event e(__func__);
  
  assertion(_hasComputedMapping);
  assertion(input()->getDimensions() == output()->getDimensions(),
            input()->getDimensions(), output()->getDimensions());

  PetscErrorCode ierr = 0;
  KSPConvergedReason convReason;
  auto& inValues = input()->data(inputDataID)->values();
  auto& outValues = output()->data(outputDataID)->values();

  int valueDim = input()->data(inputDataID)->getDimensions();
  assertion(valueDim == output()->data(outputDataID)->getDimensions(),
            valueDim, output()->data(outputDataID)->getDimensions());

  int localPolyparams = utils::Parallel::getProcessRank() > 0 ? 0 : polyparams; // Set localPolyparams only when root rank

  if (getConstraint() == CONSERVATIVE) {
    petsc::Vector au(_matrixA, "au", petsc::Vector::RIGHT);
    petsc::Vector in(_matrixA, "in");
    
    // Fill input from input data values
    for (int dim = 0; dim < valueDim; dim++) {
      printMappingInfo(inputDataID, dim);
  
      for (size_t i = 0; i < input()->vertices().size(); i++ ) {
        int globalIndex = input()->vertices()[i].getGlobalIndex();
        VecSetValue(in, globalIndex, inValues[i*valueDim + dim], INSERT_VALUES); // Dies besser als VecSetValuesLocal machen
      }
      in.assemble();

      // Gets the petsc::vector for the given combination of outputData, inputData and dimension
      // If none created yet, create one, based on _matrixC
      petsc::Vector& out = std::get<0>(
        previousSolution.emplace(std::piecewise_construct,
                                 std::forward_as_tuple(inputDataID + outputDataID * 10 + dim * 100),
                                 std::forward_as_tuple(_matrixC, "out"))
        )->second;
      
      ierr = MatMultTranspose(_matrixA, in, au); CHKERRV(ierr);
      utils::Event eSolve("PetRBF.solve.conservative");
      ierr = KSPSolve(_solver, au, out); CHKERRV(ierr);
      eSolve.stop();
      
      PetscInt iterations;
      KSPGetIterationNumber(_solver, &iterations);
      utils::EventRegistry::setProp("PetRBF.its.conservative", iterations);
      
      ierr = KSPGetConvergedReason(_solver, &convReason); CHKERRV(ierr);
      if (convReason < 0) {
        KSPView(_solver, PETSC_VIEWER_STDOUT_WORLD);
        ERROR("RBF linear system has not converged.");
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
    petsc::Vector a(_matrixQ, "a", petsc::Vector::RIGHT); // holds the solution of the LS polynom
        
    ierr = VecSetLocalToGlobalMapping(in, _ISmapping); CHKERRV(ierr);
    const PetscScalar *vecArray;

    // For every data dimension, perform mapping
    for (int dim=0; dim < valueDim; dim++) {
      printMappingInfo(inputDataID, dim);
      
      // Fill input from input data values
      int count = 0;
      for (const auto& vertex : input()->vertices()) {
        ierr = VecSetValueLocal(in, vertex.getGlobalIndex()+polyparams, inValues[count*valueDim + dim], INSERT_VALUES); CHKERRV(ierr); // evtl. besser als VecSetValuesLocal
        count++;
      }
      in.assemble();
            
      if (_polynomial == Polynomial::SEPARATE) {
        KSPSolve(_QRsolver, in, a);
        VecScale(a, -1);
        MatMultAdd(_matrixQ, a, in, in); // Subtract the polynomial from the input values
      }
      
      petsc::Vector& p = std::get<0>(  // Save and reuse the solution from the previous iteration
        previousSolution.emplace(std::piecewise_construct,
                                 std::forward_as_tuple(inputDataID + outputDataID * 10 + dim * 100),
                                 std::forward_as_tuple(_matrixC, "p"))
        )->second;

      utils::Event eSolve("PetRBF.solve.consistent");
      ierr = KSPSolve(_solver, in, p); CHKERRV(ierr);
      eSolve.stop();
            
      PetscInt iterations;
      KSPGetIterationNumber(_solver, &iterations);
      utils::EventRegistry::setProp("PetRBF.its.consistent", iterations);

      ierr = KSPGetConvergedReason(_solver, &convReason); CHKERRV(ierr);
      if (convReason < 0) {
        KSPView(_solver, PETSC_VIEWER_STDOUT_WORLD);
        ERROR("RBF linear system has not converged.");
      }

      ierr = MatMult(_matrixA, p, out); CHKERRV(ierr);

      if (useRescaling and (_polynomial == Polynomial::SEPARATE)) {
        petsc::Vector temp(_matrixA);
        ierr = MatMult(_matrixA, rescalingCoeffs, temp); CHKERRV(ierr);
        ierr = VecPointwiseDivide(out, out, temp); CHKERRV(ierr);
      }
      
      if (_polynomial == Polynomial::SEPARATE) {
        ierr = VecScale(a, -1); // scale it back, so wie add the polynom
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
void PetRadialBasisFctMapping<RADIAL_BASIS_FUNCTION_T>::tagMeshFirstRound()
{
  TRACE();
  mesh::PtrMesh filterMesh, otherMesh;
  if (getConstraint() == CONSISTENT){
    filterMesh = input();
    otherMesh = output();
  }
  else {
    assertion(getConstraint() == CONSERVATIVE, getConstraint());
    filterMesh = output();
    otherMesh = input();
  }

  for(mesh::Vertex& v : filterMesh->vertices()){
    bool isInside = true;
    if(otherMesh->vertices().size()==0) isInside = false; //ranks not at the interface should never hold interface vertices
    if(_basisFunction.hasCompactSupport()){
      for (int d=0; d<getDimensions(); d++) {
        if (v.getCoords()[d] < otherMesh->getBoundingBox()[d].first - _basisFunction.getSupportRadius() ||
            v.getCoords()[d] > otherMesh->getBoundingBox()[d].second + _basisFunction.getSupportRadius() ) {
          isInside = false;
        }
      }
    }
    if(isInside) v.tag();
  }
}

template<typename RADIAL_BASIS_FUNCTION_T>
void PetRadialBasisFctMapping<RADIAL_BASIS_FUNCTION_T>::tagMeshSecondRound()
{
  TRACE();

  if(not _basisFunction.hasCompactSupport()) return; //tags should not be changed

  mesh::Mesh::BoundingBox bb(getDimensions(), std::make_pair(std::numeric_limits<double>::max(), std::numeric_limits<double>::lowest()));

  mesh::PtrMesh mesh; //the mesh we want to filter

  if (getConstraint() == CONSISTENT){
    mesh = input();
  }
  else {
    assertion(getConstraint() == CONSERVATIVE, getConstraint());
    mesh = output();
  }

  // construct bounding box around all owned vertices
  for(mesh::Vertex& v : mesh->vertices()){
    if(v.isOwner()){
      assertion(v.isTagged());
      for (int d=0; d<getDimensions(); d++) {
        if(v.getCoords()[d] < bb[d].first) bb[d].first = v.getCoords()[d];
        if(v.getCoords()[d] > bb[d].second) bb[d].second = v.getCoords()[d];
      }
    }
  }
  // tag according to bounding box
  for(mesh::Vertex& v : mesh->vertices()){
    bool isInside = true;
    for (int d=0; d<getDimensions(); d++) {
      if (v.getCoords()[d] < bb[d].first - _basisFunction.getSupportRadius() ||
          v.getCoords()[d] > bb[d].second + _basisFunction.getSupportRadius() ) {
        isInside = false;
      }
    }
    if(isInside) v.tag();
  }
}


template <typename RADIAL_BASIS_FUNCTION_T>
void PetRadialBasisFctMapping<RADIAL_BASIS_FUNCTION_T>::printMappingInfo(int inputDataID, int dim) const
{
  const std::string constraintName = getConstraint() == CONSERVATIVE ? "conservative" : "consistent";
  const std::string polynomialName = _polynomial == Polynomial::ON ? "on" : _polynomial == Polynomial::OFF ? "off" : "separate";

  INFO("Mapping " << input()->data(inputDataID)->getName() << " " << constraintName
                  << " from " << input()->getName() << " (ID " << input()->getID() << ")"
                  << " to " << output()->getName() << " (ID " << output()->getID() << ") "
                  << "for dimension " << dim << ") with polynomial set to " << polynomialName);
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
