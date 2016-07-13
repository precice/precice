#include "Petsc.hpp"
#include "utils/Globals.hpp"

#ifndef PRECICE_NO_PETSC
#include "petsc.h"
#endif // not PRECICE_NO_PETSC

namespace precice {
namespace utils {

tarch::logging::Log Petsc:: _log ( "precice::utils::Petsc" );

bool Petsc::weInitialized = false;

void Petsc::initialize
(
  int*               argc,
  char***            argv )
{
  preciceTrace ("initialize()");
#ifndef PRECICE_NO_PETSC
  PetscBool petscIsInitialized;
  PetscInitialized(&petscIsInitialized);
  if (not petscIsInitialized) {
    PetscErrorCode ierr;
    ierr = PetscInitialize(argc, argv, "", nullptr); CHKERRV(ierr);
    weInitialized = true;
    PetscPushErrorHandler(&PetscMPIAbortErrorHandler, nullptr);
  }
#endif // not PRECICE_NO_PETSC
}

void Petsc::finalize()
{
#ifndef PRECICE_NO_PETSC
  PetscBool petscIsInitialized;
  PetscInitialized(&petscIsInitialized);
  if (petscIsInitialized and weInitialized) {
    PetscFinalize();
  }
#endif // not PRECICE_NO_PETSC
}
}} // namespace precice, utils


#ifndef PRECICE_NO_PETSC

#include <string>
#include <limits>
#include <random>
#include "petscviewer.h"
#include "petscdraw.h"

namespace precice {
namespace utils {
namespace petsc {

Vector::Vector(MPI_Comm comm, std::string name)
{
  PetscErrorCode ierr = 0;
  ierr = VecCreate(comm, &vector); CHKERRV(ierr);
  setName(name);
}

Vector::Vector(Vec &v, std::string name)
{
  VecCopy(v, vector);
  setName(name);
}

Vector::Vector(Vector &v, std::string name)
{
  PetscErrorCode ierr = 0;
  ierr = VecDuplicate(v.vector, &vector); CHKERRV(ierr);
  setName(name);
}

Vector::Vector(Mat &m, std::string name, LEFTRIGHT type)
{
  // MatGetVecs is deprecated, we keep it due to the old PETSc version at the SuperMUC.
  PetscErrorCode ierr = 0;
  if (type == LEFTRIGHT::LEFT) {
    ierr = MatGetVecs(m, nullptr, &vector); CHKERRV(ierr); // a vector with the same number of rows
  }
  else {
    ierr = MatGetVecs(m, &vector, nullptr); CHKERRV(ierr); // a vector with the same number of cols
  }
  setName(name);
}

Vector::Vector(Matrix &m, std::string name, LEFTRIGHT type) :
  Vector(m.matrix, name, type)
{}

Vector::~Vector()
{
  PetscErrorCode ierr = 0;
  PetscBool petscIsInitialized;
  PetscInitialized(&petscIsInitialized);
  if (petscIsInitialized) // If PetscFinalize is called before ~Vector
    ierr = VecDestroy(&vector); CHKERRV(ierr);
}

Vector::operator Vec&()
{
  return vector;
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

int Vector::getLocalSize()
{
  PetscInt size;
  VecGetLocalSize(vector, &size);
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
    a[i - range_start] = (i + start) * step_size;
  }
  VecRestoreArray(vector, &a);
}

void Vector::fill_with_randoms()
{
  PetscErrorCode ierr = 0;
  PetscRandom rctx;

  std::random_device rd;
  std::uniform_real_distribution<double> dist(0, 1);

  PetscRandomCreate(PETSC_COMM_WORLD, &rctx);
  PetscRandomSetType(rctx, PETSCRAND48);
  PetscRandomSetSeed(rctx, dist(rd));
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


std::pair<PetscInt, PetscInt> Vector::ownerRange()
{
  PetscInt range_start, range_end;
  VecGetOwnershipRange(vector, &range_start, &range_end);
  return std::make_pair(range_start, range_end);
}
  
void Vector::write(std::string filename)
{
  PetscErrorCode ierr = 0;
  PetscViewer fd;
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename.c_str(), FILE_MODE_WRITE, &fd); CHKERRV(ierr);
  VecView(vector, fd); CHKERRV(ierr);
  PetscViewerDestroy(&fd);
}

void Vector::read(std::string filename)
{
   PetscErrorCode ierr = 0;
   PetscViewer fd;
   PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename.c_str(), FILE_MODE_READ, &fd); CHKERRV(ierr);
   VecLoad(vector, fd); CHKERRV(ierr); CHKERRV(ierr);
   PetscViewerDestroy(&fd);
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
  PetscBool petscIsInitialized;
  PetscInitialized(&petscIsInitialized);
  if (petscIsInitialized) // If PetscFinalize is called before ~Matrix
    ierr = MatDestroy(&matrix); CHKERRV(ierr);
}

Matrix::operator Mat&()
{
  return matrix;
}

void Matrix::assemble(MatAssemblyType type)
{
  PetscErrorCode ierr = 0;
  ierr = MatAssemblyBegin(matrix, type); CHKERRV(ierr);
  ierr = MatAssemblyEnd(matrix, type); CHKERRV(ierr);
}

void Matrix::init(PetscInt localRows, PetscInt localCols, PetscInt globalRows, PetscInt globalCols,
                  MatType type, bool doSetup)
{
  PetscErrorCode ierr = 0;
  if (type != nullptr) {
    ierr = MatSetType(matrix, type); CHKERRV(ierr);
  }
  ierr = MatSetSizes(matrix, localRows, localCols, globalRows, globalCols); CHKERRV(ierr);
  ierr = MatSetFromOptions(matrix); CHKERRV(ierr);
  if (doSetup)
    ierr = MatSetUp(matrix); CHKERRV(ierr);
}

void Matrix::reset()
{
  PetscErrorCode ierr = 0;
  std::string name = getName();
  ierr = MatDestroy(&matrix); CHKERRV(ierr);
  ierr = MatCreate(communicator, &matrix); CHKERRV(ierr);
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

MatInfo Matrix::getInfo(MatInfoType flag)
{
  MatInfo info;
  MatGetInfo(matrix, flag, &info);
  return info;
}

void Matrix::setValue(PetscInt row, PetscInt col, PetscScalar value)
{
  PetscErrorCode ierr = 0;
  ierr = MatSetValue(matrix, row, col, value, INSERT_VALUES); CHKERRV(ierr);
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

std::pair<PetscInt, PetscInt> Matrix::getSize()
{
  PetscInt m, n;
  MatGetSize(matrix, &m, &n);
  return std::make_pair(m, n);
}

std::pair<PetscInt, PetscInt> Matrix::ownerRange()
{
  PetscInt range_start, range_end;
  MatGetOwnershipRange(matrix, &range_start, &range_end);
  return std::make_pair(range_start, range_end);
}

std::pair<PetscInt, PetscInt> Matrix::ownerRangeColumn()
{
  PetscInt range_start, range_end;
  MatGetOwnershipRangeColumn(matrix, &range_start, &range_end);
  return std::make_pair(range_start, range_end);
}

void Matrix::write(std::string filename)
{
  PetscErrorCode ierr = 0;
  PetscViewer fd;
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename.c_str(), FILE_MODE_WRITE, &fd); CHKERRV(ierr);
  MatView(matrix, fd); CHKERRV(ierr);
  PetscViewerDestroy(&fd);
}

void Matrix::read(std::string filename)
{
   PetscErrorCode ierr = 0;
   PetscViewer fd;
   PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename.c_str(), FILE_MODE_READ, &fd); CHKERRV(ierr);
   MatLoad(matrix, fd); CHKERRV(ierr); CHKERRV(ierr);
   PetscViewerDestroy(&fd);
}

void Matrix::view()
{
  PetscErrorCode ierr = 0;
  PetscViewer viewer;
  ierr = PetscViewerCreate(communicator, &viewer); CHKERRV(ierr);
  ierr = PetscViewerSetType(viewer, PETSCVIEWERASCII); CHKERRV(ierr); 
  ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_DENSE); CHKERRV(ierr);
  ierr = MatView(matrix, viewer); CHKERRV(ierr);
  ierr = PetscViewerDestroy(&viewer); CHKERRV(ierr); 
}

void Matrix::viewDraw()
{
  PetscErrorCode ierr = 0;
  PetscViewer viewer;
  PetscDraw draw;
  ierr = PetscViewerCreate(communicator, &viewer); CHKERRV(ierr);
  ierr = PetscViewerSetType(viewer, PETSCVIEWERDRAW); CHKERRV(ierr); 
  ierr = MatView(matrix, viewer); CHKERRV(ierr);
  ierr = PetscViewerDrawGetDraw(viewer, 0, &draw); CHKERRV(ierr);
  ierr = PetscDrawSetPause(draw, -2); CHKERRV(ierr); // pause on destroy
  ierr = PetscViewerDestroy(&viewer); CHKERRV(ierr);
  ierr = PetscDrawDestroy(&draw); CHKERRV(ierr); 
}

}}} // namespace precice, utils, petsc

#endif // PRECICE_NO_PETSC




