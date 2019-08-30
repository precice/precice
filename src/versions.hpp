#pragma once

#define PRECICE_VERSION "1.3.0"

#ifndef PRECICE_NO_PETSC
  #define PETSC_MAJOR 0
  #define PETSC_MINOR 0
#endif

namespace precice {
    extern const char * preciceRevision;
    extern const char * versionInformation;
}
