#ifndef PRECICE_NO_PETSC

#include <string>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <limits>
#include "petnum.hpp"

#include "petscviewer.h" 

namespace petsc {


Vector::Vector(MPI_Comm comm, std::string name)
{
  PetscErrorCode ierr = 0;
  ierr = VecCreate(PETSC_COMM_WORLD, &vector); CHKERRV(ierr);
  setName(name);
}

Vector::Vector(Vec &v, std::string name)
{
  vector = v;
  setName(name);
}

Vector::Vector(Vector &v, std::string name)
{
  PetscErrorCode ierr = 0;
  ierr = VecDuplicate(v.vector, &vector); CHKERRV(ierr);
  setName(name);
}

Vector::Vector(Mat &m, std::string name)
{
  PetscErrorCode ierr = 0;
  ierr = MatGetVecs(m, &vector, nullptr); CHKERRV(ierr); // a vector with the same number of rows
  setName(name);
}

Vector::Vector(Matrix &m, std::string name)
{
  PetscErrorCode ierr = 0;
  ierr = MatGetVecs(m.matrix, nullptr, &vector); CHKERRV(ierr); // a vector with the same number of rows
  setName(name);
}

// Vector::Vector(const Vector &v, std::string name)
// {
//   PetscErrorCode ierr = 0;
//   ierr = VecDuplicate(v.vector, &vector); CHKERRV(ierr);
//   if (!name.empty())
//     ierr = PetscObjectSetName( (PetscObject)vector, name.c_str()); CHKERRV(ierr);    
// }

Vector::~Vector()
{
  PetscErrorCode ierr = 0;
  ierr = VecDestroy(&vector); CHKERRV(ierr);
}

void Vector::init(PetscInt rows)
{
  PetscErrorCode ierr = 0;
  ierr = VecSetSizes(vector, PETSC_DECIDE, rows); CHKERRV(ierr);
  ierr = VecSetFromOptions(vector); CHKERRV(ierr);
}

void Vector::setName(std::string name)
{
  PetscErrorCode ierr = 0;
  ierr = PetscObjectSetName((PetscObject) vector, name.c_str()); CHKERRV(ierr); 
}

std::string Vector::getName()
{
  const char *cstr;
  PetscObjectGetName((PetscObject) vector, &cstr); 
  return cstr;    
}


int Vector::getSize()
{
  PetscInt size;
  VecGetSize(vector, &size);
  return size;
}

void Vector::setValue(PetscInt row, PetscScalar value)
{
  PetscErrorCode ierr = 0;
  ierr = VecSetValue(vector, row, value, INSERT_VALUES); CHKERRV(ierr);
}

void Vector::arange(double start, double stop)
{
  PetscErrorCode ierr = 0;
  PetscScalar *a;
  PetscInt range_start, range_end, size;
  VecGetSize(vector, &size);
  VecGetOwnershipRange(vector, &range_start, &range_end);
  double step_size = (stop-start) / size;
  ierr = VecGetArray(vector, &a); CHKERRV(ierr); 
  for (PetscInt i = range_start; i < range_end; i++) {
    // a[i - range_start] = i + start;
    a[i - range_start] = (i + start) * step_size;
  }
  VecRestoreArray(vector, &a);
}

void Vector::fill_with_randoms()
{
  PetscErrorCode ierr = 0;
  PetscRandom rctx;

  unsigned long seed = (double) std::rand()/RAND_MAX  * std::numeric_limits<unsigned long>::max();
    
  PetscRandomCreate(PETSC_COMM_WORLD, &rctx);
  PetscRandomSetType(rctx, PETSCRAND48);
  PetscRandomSetSeed(rctx, seed);
  PetscRandomSeed(rctx);     
  ierr = VecSetRandom(vector, rctx); CHKERRV(ierr);
  PetscRandomDestroy(&rctx);
}

void Vector::sort() 
{
  // will not work on multiple processors as expected since only local partion will be sorted.
  PetscErrorCode ierr = 0;
  PetscInt size;
  PetscReal *a;
  ierr = VecGetArray(vector, &a); CHKERRV(ierr);
  ierr = VecGetSize(vector, &size);
  ierr = PetscSortReal(size, a); CHKERRV(ierr);
  ierr = VecRestoreArray(vector, &a); CHKERRV(ierr);
}

void Vector::assemble()
{
  PetscErrorCode ierr = 0;
  ierr = VecAssemblyBegin(vector); CHKERRV(ierr); 
  ierr = VecAssemblyEnd(vector); CHKERRV(ierr); 
}
  
void Vector::write(std::string filename)
{
  PetscViewer viewer;
  PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename.c_str(), &viewer);
  PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
  VecView(vector, viewer);
  PetscViewerDestroy(&viewer);
        
}

void Vector::view()
{
  PetscErrorCode ierr;
  ierr = VecView(vector, PETSC_VIEWER_STDOUT_WORLD); CHKERRV(ierr);
}


/////////////////////////////////////////////////////////////////////////

Matrix::Matrix(MPI_Comm comm, std::string name)
  :
  communicator(comm)
{
  PetscErrorCode ierr = 0;
  ierr = MatCreate(communicator, &matrix); CHKERRV(ierr);
  setName(name);
}
  

Matrix::~Matrix()
{
  PetscErrorCode ierr = 0;
  ierr = MatDestroy(&matrix); CHKERRV(ierr);
}

void Matrix::reset()
{
  PetscErrorCode ierr = 0;
  std::string name = getName();
  ierr = MatDestroy(&matrix); CHKERRV(ierr);
  ierr = MatCreate(communicator, &matrix); CHKERRV(ierr);
  setName(name);
}

void Matrix::create(size_t rows, size_t cols, std::string name)
{
  PetscErrorCode ierr = 0;
  ierr = MatSetType(matrix, MATDENSE);
  ierr = MatSetSizes(matrix, PETSC_DECIDE, PETSC_DECIDE, rows, cols); CHKERRV(ierr);
  ierr = MatSetFromOptions(matrix); CHKERRV(ierr);
  ierr = MatSetUp(matrix); CHKERRV(ierr); 
  // ierr = MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY); CHKERRV(ierr); 
  // ierr = MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY); CHKERRV(ierr); 
  setName(name);
}

void Matrix::setName(std::string name)
{
  PetscErrorCode ierr = 0;
  ierr = PetscObjectSetName((PetscObject) matrix, name.c_str()); CHKERRV(ierr); 
}

std::string Matrix::getName()
{
  const char *cstr;
  PetscObjectGetName((PetscObject) matrix, &cstr); 
  return cstr;    
}


void Matrix::assemble(MatAssemblyType type)
{
  PetscErrorCode ierr = 0;
  ierr = MatAssemblyBegin(matrix, type); CHKERRV(ierr);
  ierr = MatAssemblyEnd(matrix, type); CHKERRV(ierr);
}


void Matrix::fill_with_randoms()
{
  PetscErrorCode ierr = 0;
  PetscRandom rctx;

  unsigned long seed = (double) std::rand()/RAND_MAX * std::numeric_limits<unsigned long>::max();
    
  PetscRandomCreate(communicator, &rctx);
  PetscRandomSetType(rctx, PETSCRAND48);
  PetscRandomSetSeed(rctx, seed);
  PetscRandomSeed(rctx);     
  ierr = MatSetRandom(matrix, rctx); CHKERRV(ierr);
  PetscRandomDestroy(&rctx);
}

void Matrix::set_column(Vector &v, int col)
{
  PetscErrorCode ierr = 0;
  const PetscScalar *vec;
  PetscInt range_start, range_end;
  VecGetOwnershipRange(v.vector, &range_start, &range_end);
  PetscInt irow[range_end - range_start];
  for (PetscInt i = range_start; i < range_end; i++) {
    irow[i - range_start] = i;
  }
      
  ierr = VecGetArrayRead(v.vector, &vec); CHKERRV(ierr);
  ierr = MatSetValues(matrix, range_end - range_start, irow, 1, &col, vec, INSERT_VALUES); CHKERRV(ierr);
  ierr = VecRestoreArrayRead(v.vector, &vec); CHKERRV(ierr);
  ierr = MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY); CHKERRV(ierr); 
  ierr = MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY); CHKERRV(ierr); 
}

  
void Matrix::write(std::string filename)
{
  PetscViewer viewer;
  std::remove(filename.c_str()); // ugly
  PetscViewerASCIIOpen(communicator, filename.c_str(), &viewer);
  PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
  MatView(matrix, viewer);
  PetscViewerDestroy(&viewer);
        
}

void Matrix::view()
{
  PetscErrorCode ierr;
  PetscViewer viewer;
  ierr = PetscViewerCreate(communicator, &viewer); CHKERRV(ierr);
  ierr = PetscViewerSetType(viewer, PETSCVIEWERASCII); CHKERRV(ierr); 
  ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_DENSE); CHKERRV(ierr);
  ierr = MatView(matrix, viewer); CHKERRV(ierr);
  ierr = PetscViewerDestroy(&viewer); CHKERRV(ierr); 
}

}

#endif
