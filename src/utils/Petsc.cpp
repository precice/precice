#include "Petsc.hpp"
#include "utils/Globals.hpp"

#ifndef PRECICE_NO_PETSC
#include "petsc.h"
#endif // not PRECICE_NO_PETSC

namespace precice {
namespace utils {

tarch::logging::Log Petsc:: _log ( "precice::utils::Petsc" );


void Petsc:: initialize
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
  }
#endif // not PRECICE_NO_PETSC
}

void Petsc:: finalize()
{
#ifndef PRECICE_NO_PETSC
  PetscBool petscIsInitialized;
  PetscInitialized(&petscIsInitialized);
  if (petscIsInitialized) {
    PetscFinalize();
  }
#endif // not PRECICE_NO_PETSC
}

}} // precice, utils


