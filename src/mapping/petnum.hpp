#pragma once

#include <string>

#include "petscvec.h"
#include "petscmat.h"

namespace petsc {

class Vector
{
public:
  Vec vector;
  
  Vector(size_t size, std::string name = "");
  Vector(Vec &v, std::string name = "");
  Vector(const Vector &v, std::string name = "");

  ~Vector();

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
  Matrix(MPI_Comm comm);

  ~Matrix();

  void create(size_t rows, size_t cols, std::string name = "");

  void setName(std::string name);
  
  void assemble();
  
  void fill_with_randoms();
  
  void set_column(Vector &v, int col);

  void write(std::string filename);
  
  void view();

};
}
