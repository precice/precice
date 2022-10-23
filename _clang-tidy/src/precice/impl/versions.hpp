#pragma once

#define PRECICE_VERSION "2.4.0"

#ifndef PRECICE_NO_PETSC
#define PETSC_MAJOR 3
#define PETSC_MINOR 12
#endif

namespace precice {
extern char const *const preciceRevision;
extern char const *const versionInformation;
} // namespace precice
