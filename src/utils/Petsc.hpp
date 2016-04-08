#pragma once

#include "logging/Logger.hpp"

namespace precice {
namespace utils {

/// Utility class for managing PETSc operations.
class Petsc
{
public:
  /**
   * @brief Initializes the Petsc environment.
   *
   * @param[in] argc Parameter count, passed to PetscInitialize
   * @param[in] argv Parameter values, passed to PetscInitialize
   */
  static void initialize (
    int*               argc,
    char***            argv);

  /// Finalizes Petsc environment.
  static void finalize();

private:
  
  /// Whether we have initialized Petsc or if it was initialized by an application calling us.
  static bool weInitialized;
  
  static logging::Logger _log;
};
}} // namespace precice, utils


#ifndef PRECICE_NO_PETSC

#include <string>
#include <utility>

#include "petscvec.h"
#include "petscmat.h"

namespace precice {
namespace utils {
namespace petsc {

class Matrix;

class Vector
{
public:
  Vec vector;

  enum LEFTRIGHT { LEFT, RIGHT };
  
  /// Creates a new vector on the given MPI communicator.
  Vector(MPI_Comm comm = PETSC_COMM_WORLD, std::string name = "");

  /// Use Vec v as vector.
  Vector(Vec &v, std::string name = "");

  /// Duplicates type, row layout etc. (not values) of v.
  Vector(Vector &v, std::string name = "");  

  /// Constructs a vector with the same number of rows (default) or columns.
  Vector(Mat &m, std::string name = "", LEFTRIGHT type = LEFT);

  /// Constructs a vector with the same number of rows (default) or columns.
  Vector(Matrix &m, std::string name = "", LEFTRIGHT type = LEFT);

  /// Delete copy and assignement constructor
  /** Copying and assignement of this class would involve copying the pointer to
   * the PETSc object and finallly cause double destruction of it.
   */
  Vector(const Vector&) = delete;
  Vector& operator=(const Vector&) = delete;

  ~Vector();

  /// Sets the size and calls VecSetFromOptions
  void init(PetscInt rows);

  void setName(std::string name);
  std::string getName();

  int getSize();

  int getLocalSize();
  
  void setValue(PetscInt row, PetscScalar value);

  void arange(double start, double stop);

  void fill_with_randoms();

  void sort();

  void assemble();

  /// Returns a pair that mark the beginning and end of the vectors ownership range. Use first und second to access.
  std::pair<PetscInt, PetscInt> ownerRange();

  /// Writes the vector to the PETSc binary format
  void write(std::string filename);

  /// Reads the vector from the PETSc binary format
  void read(std::string filename);

  void view();
};

  
class Matrix
{
public:
  Mat matrix;

  MPI_Comm communicator;

  /// Delete copy and assignement constructor
  /* Copying and assignement of this class would involve copying the pointer to
   * the PETSc object and finallly cause double destruction of it.
   */
  Matrix(const Matrix&) = delete;
  Matrix& operator=(const Matrix&) = delete;

  Matrix(MPI_Comm comm = PETSC_COMM_WORLD, std::string name = "");

  ~Matrix();

  void assemble(MatAssemblyType type = MAT_FINAL_ASSEMBLY);
    
  /// Initializes matrix of given size and type
  void init(PetscInt localRows, PetscInt localCols, PetscInt globalRows, PetscInt globalCols, MatType type = nullptr);

  /// Destroys and recreates the matrix on the same communicator
  void reset();
  
  void setName(std::string name);
  std::string getName();

  /// Get the MatInfo struct for the matrix.
  /* See http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Mat/MatInfo.html for description of fields.
   */
  MatInfo getInfo(MatInfoType flag);
  
  void setValue(PetscInt row, PetscInt col, PetscScalar value);
  
  void fill_with_randoms();
  
  void set_column(Vector &v, int col);

  std::pair<PetscInt, PetscInt> getSize();
  
  /// Returns a pair that mark the beginning and end of the matrix' ownership range.
  std::pair<PetscInt, PetscInt> ownerRange();
  
  /// Returns a pair that mark the beginning and end of the matrix' column ownership range.
  std::pair<PetscInt, PetscInt> ownerRangeColumn();
  
  /// Writes the matrix to PETSc the binary format
  void write(std::string filename);

  /// Reads the matrix from PETSc the binary format
  void read(std::string filename);

  /// Prints the matrix
  void view();

  void viewDraw();
  
};

}}} // namespace precice, utils, petsc

#endif // PRECICE_NO_PETSC




