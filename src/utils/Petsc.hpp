#pragma once

#include "logging/Logger.hpp"
#include "utils/Parallel.hpp"

namespace precice {
namespace logging {
class Logger;
} // namespace logging

namespace utils {

/// Utility class for managing PETSc operations.
class Petsc {
public:
  /**
   * @brief Initializes the Petsc environment.
   *
   * @param[in] argc Parameter count, passed to PetscInitialize
   * @param[in] argv Parameter values, passed to PetscInitialize
   * @param[in] comm The communicator to Initialize PETSc on
   */
  static void initialize(
      int *                         argc,
      char ***                      argv,
      utils::Parallel::Communicator comm);

  /// Finalizes Petsc environment.
  static void finalize();

private:
  /// Whether we have initialized Petsc or if it was initialized by an application calling us.
  static bool weInitialized;

  static logging::Logger _log;
};
} // namespace utils
} // namespace precice

#ifndef PRECICE_NO_PETSC

#include <string>
#include <utility>
#include "petscao.h"
#include "petscis.h"
#include "petscksp.h"
#include "petscmat.h"
#include "petscvec.h"

namespace precice {
namespace utils {
/// PETSc related utilities
namespace petsc {

enum VIEWERFORMAT { ASCII,
                    BINARY };

class Matrix;

class Vector {
public:
  Vec vector = nullptr;

  enum LEFTRIGHT { LEFT,
                   RIGHT };

  /// Creates a new vector on the given MPI communicator.
  explicit Vector(const std::string &name = "");

  /** Copy construction from another vector
   * Duplicates the vector and copies the name
   */
  Vector(const Vector &other);

  /** Copy assignement
   * Destroys the current vector and takes ownership of the other.
   */
  Vector &operator=(const Vector &other);

  /** Move construction
   * Takes ownership of the other vector.
   */
  Vector(Vector &&other);

  /** Move assignement
   * Destroys the current vector and takes ownership of the other.
   */
  Vector &operator=(Vector &&other);

  /** Constructs the object from another Vec
   * Takes ownership of the other Vec
   */
  Vector(Vec &other, const std::string &name = "");

  ~Vector();

  ///@name Allocation
  ///@{

  /// Allocates a new vector on the given MPI communicator.
  static Vector allocate(const std::string &name = "");

  /** Allocated an uninitialized vector of identical shape.
   * Duplicates type, row layout etc. (not values) of v.
   */
  static Vector allocate(Vector &other, const std::string &name = "");

  /// Allocated an uninitialized vector of identical shape.
  static Vector allocate(Vec &other, const std::string &name = "");

  /// Allocates a vector with the same number of rows (default) or columns.
  static Vector allocate(Matrix &m, const std::string &name = "", LEFTRIGHT type = LEFT);

  /// Allocates a vector with the same number of rows (default) or columns.
  static Vector allocate(Mat &m, const std::string &name = "", LEFTRIGHT type = LEFT);

  ///@}

  /// Swaps the ownership of two vectors
  void swap(Vector &other) noexcept;

  /// Enables implicit conversion into a reference to a PETSc Vec type
  operator Vec &();

  /// Sets the size and calls VecSetFromOptions
  void init(PetscInt rows);

  PetscInt getSize() const;

  PetscInt getLocalSize() const;

  void setValue(PetscInt row, PetscScalar value);

  void arange(double start, double stop);

  void fillWithRandoms();

  /// Sorts the LOCAL partion of the vector
  void sort();

  void assemble();

  /// Returns a pair that mark the beginning and end of the vectors ownership range. Use first und second to access.
  std::pair<PetscInt, PetscInt> ownerRange() const;

  /// Writes the vector to file.
  void write(std::string filename, VIEWERFORMAT format = ASCII) const;

  /// Reads the vector from file.
  void read(std::string filename, VIEWERFORMAT format = ASCII);

  /// Prints the vector
  void view() const;
};

void swap(Vector &lhs, Vector &rhs) noexcept;

class Matrix {
public:
  Mat matrix = nullptr;

  /// Delete copy and assignment constructor
  /** Copying and assignment of this class would involve copying the pointer to
      the PETSc object and finallly cause double destruction of it.
  */
  Matrix(const Matrix &) = delete;
  Matrix &operator=(const Matrix &) = delete;

  explicit Matrix(std::string name = "");

  /// Move constructor, use the implicitly declared.
  Matrix(Matrix &&) = default;
  Matrix &operator=(Matrix &&) = default;

  ~Matrix();

  /// Enables implicit conversion into a reference to a PETSc Mat type
  operator Mat &();

  void assemble(MatAssemblyType type = MAT_FINAL_ASSEMBLY);

  /// Initializes matrix of given size and type
  /** @param[in] localRows,localCols The number of rows/cols that are local to the processor
      @param[in] globalRows,globalCols The number of global rows/cols.
      @param[in] type PETSc type of the matrix
      @param[in] doSetup Call MatSetup(). Not calling MatSetup can have performance gains when using preallocation
  */
  void init(PetscInt localRows, PetscInt localCols, PetscInt globalRows, PetscInt globalCols,
            MatType type = nullptr, bool doSetup = true);

  /// Destroys and recreates the matrix on the same communicator
  void reset();

  /// Get the MatInfo struct for the matrix.
  /** See http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Mat/MatInfo.html for description of fields. */
  MatInfo getInfo(MatInfoType flag) const;

  void setValue(PetscInt row, PetscInt col, PetscScalar value);

  void fillWithRandoms();

  void setColumn(Vector &v, PetscInt col);

  /// Returns (rows, cols) global size
  std::pair<PetscInt, PetscInt> getSize() const;

  /// Returns (rows, cols) local size
  std::pair<PetscInt, PetscInt> getLocalSize() const;

  /// Returns a pair that mark the beginning and end of the matrix' ownership range.
  std::pair<PetscInt, PetscInt> ownerRange() const;

  /// Returns a pair that mark the beginning and end of the matrix' column ownership range.
  std::pair<PetscInt, PetscInt> ownerRangeColumn() const;

  /// Returns the block size of the matrix
  PetscInt blockSize() const;

  /// Writes the matrix to file.
  void write(std::string filename, VIEWERFORMAT format = ASCII) const;

  /// Reads the matrix from file, stored in PETSc binary format
  void read(std::string filename);

  /// Prints the matrix
  void view() const;

  /// Graphically draws the matrix structure
  void viewDraw() const;
};

class KSPSolver {
public:
  KSP ksp;

  /// Delete copy and assignment constructor
  /** Copying and assignment of this class would involve copying the pointer to
      the PETSc object and finallly cause double destruction of it.
  */
  KSPSolver(const KSPSolver &) = delete;
  KSPSolver &operator=(const KSPSolver &) = delete;

  explicit KSPSolver(std::string name = "");

  /// Move constructor, use the implicitly declared.
  KSPSolver(KSPSolver &&other) = default;

  ~KSPSolver();

  /// Enables implicit conversion into a reference to a PETSc KSP type
  operator KSP &();

  /// Destroys and recreates the ksp on the same communicator
  void reset();

  /// Solves the linear system, returns false it not converged
  bool solve(Vector &b, Vector &x);

  /// Solves the transposed linear system, returns false it not converged
  bool solveTranspose(Vector &b, Vector &x);

  /// Returns the iteration number of solver, either during or after the solve call.
  PetscInt getIterationNumber();
};

/// Destroys an KSP, if ksp is not null and PetscIsInitialized
void destroy(KSP *ksp);

/// Destroys an ISLocalToGlobalMapping, if IS is not null and PetscIsInitialized
void destroy(ISLocalToGlobalMapping *IS);

/// Destroys an application ordering, if ao is not null and PetscIsInitialized
void destroy(AO *ao);

} // namespace petsc
} // namespace utils
} // namespace precice

#endif // PRECICE_NO_PETSC
