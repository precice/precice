#pragma once

#include <string>

#include "petscvec.h"
#include "petscmat.h"

namespace petsc {

class Matrix;

class Vector
{
public:
  Vec vector;
  Vector(MPI_Comm comm, std::string name = "");
  Vector(Vec &v, std::string name = "");
  Vector(Vector &v, std::string name = "");

  //@brief Constructs a vector with the same number of rows.
  Vector(Mat &m, std::string name = "");

  Vector(Matrix &m, std::string name = "");

  ~Vector();

  void init(PetscInt rows);

  void setName(std::string name);
  std::string getName();

  int getSize();

  void setValue(PetscInt row, PetscScalar value);

  void arange(double start, double stop);

  void fill_with_randoms();

  void sort();

  void assemble();
  
  void write(std::string filename);

  void view();
  
  
};
  
class Matrix
{
public:
  Mat matrix;

  MPI_Comm communicator;

  Matrix() { };
  Matrix(MPI_Comm comm, std::string name = "");

  ~Matrix();

  void reset();

  void create(size_t rows, size_t cols, std::string name = "");

  void setName(std::string name);
  std::string getName();
  
  void assemble(MatAssemblyType type = MAT_FINAL_ASSEMBLY);
  
  void fill_with_randoms();
  
  void set_column(Vector &v, int col);

  void write(std::string filename);
  
  void view();

private:
  std::string _name;

};
}
