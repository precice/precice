#ifndef PRECICE_NO_PETSC
#pragma once

#include <string>
#include <utility>

#include "petscvec.h"
#include "petscmat.h"

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

  /// Constructs a vector with the same number of rows.
  Vector(Mat &m, std::string name = "");

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
    
  void view();

  void write(std::string filename);  
};

  
class Matrix
{
public:
  Mat matrix;

  MPI_Comm communicator;

  /// Delete copy and assignement constructor
  /** Copying and assignement of this class would involve copying the pointer to
   * the PETSc object and finallly cause double destruction of it.
   */
  Matrix(const Matrix&) = delete;
  Matrix& operator=(const Matrix&) = delete;

  Matrix(MPI_Comm comm = PETSC_COMM_WORLD, std::string name = "");

  ~Matrix();

  void assemble(MatAssemblyType type = MAT_FINAL_ASSEMBLY);
    
  /// Initializes matrix of given size and type
  void init(PetscInt localRows, PetscInt localCols, PetscInt globalRows, PetscInt globalCols, MatType type = nullptr);

  // Destroys and recreate the matrix on the same communicator  
  void reset();
  
  void setName(std::string name);
  std::string getName();

  void setValue(PetscInt row, PetscInt col, PetscScalar value);
  
  void fill_with_randoms();
  
  void set_column(Vector &v, int col);

  std::pair<PetscInt, PetscInt> getSize();
  
  /// Returns a pair that mark the beginning and end of the matrix' ownership range. Use first und second to access.
  std::pair<PetscInt, PetscInt> ownerRange();

  void write(std::string filename);
  
  void view();

  void viewDraw();
  
};
}

#endif
