#pragma once

// clang-format off

#define PRECICE_VERSION_MAJOR 2
#define PRECICE_VERSION_MINOR 4
#define PRECICE_VERSION_PATCH 0

#define PRECICE_VERSION "2.4.0"

// clang-format on

#define PRECICE_VERSION_GREATER_EQUAL(major, minor, patch) (             \
    (PRECICE_VERSION_MAJOR > major) ||                                   \
    (PRECICE_VERSION_MAJOR == major && PRECICE_VERSION_MINOR > minor) || \
    (PRECICE_VERSION_MAJOR == major && PRECICE_VERSION_MINOR == minor && PRECICE_VERSION_PATCH >= patch))
