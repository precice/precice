#pragma once

/**@file
 * This file contains verification macros to harden the Participant against misuse and misconfiguration.
 *
 * \attention Only include this file from ParticipantImpl.cpp
 */

//
// MESH VALIDATION
//

/** Implementation of PRECICE_VALIDATE_MESH_ID()
 *
 * @attention Do not use this macro directly!
 */
#include "utils/stacktrace.hpp"
#define PRECICE_VALIDATE_MESH_NAME_IMPL(name)                     \
  PRECICE_CHECK(_accessor->hasMesh(name),                         \
                "The mesh named \"{}\" is unknown to preCICE.{}", \
                name, _accessor->hintForMesh(name));

/** Implementation of PRECICE_REQUIRE_MESH_USE()
 *
 * @attention Do not use this macro directly!
 */
#define PRECICE_REQUIRE_MESH_USE_IMPL(name)                                                                    \
  PRECICE_VALIDATE_MESH_NAME_IMPL(name)                                                                        \
  PRECICE_CHECK(_accessor->isMeshUsed(name),                                                                   \
                "This participant does not use the mesh \"{0}\", but attempted to access it. "                 \
                "Please define a <provide-mesh name=\"{0}\" /> or <receive-mesh name=\"{0}\" from=\"...\" /> " \
                "tag in the configuration of participant \" {1}.",                                             \
                name, _accessorName);

/** Implementation of PRECICE_REQUIRE_MESH_PROVIDE()
 *
 * @attention Do not use this macro directly!
 */
#define PRECICE_REQUIRE_MESH_PROVIDE_IMPL(name)                                                \
  PRECICE_REQUIRE_MESH_USE_IMPL(name)                                                          \
  PRECICE_CHECK(_accessor->isMeshProvided(name),                                               \
                "This participant does not provide Mesh \"{0}\", but attempted to modify it. " \
                "Please add a provide-mesh tag as follows <provide-mesh name=\"{0}\" />.",     \
                name);

/** Implementation of PRECICE_REQUIRE_MESH_MODIFY()
 *
 * @attention Do not use this macro directly!
 */
#define PRECICE_REQUIRE_MESH_MODIFY_IMPL(name)                                        \
  PRECICE_REQUIRE_MESH_PROVIDE_IMPL(name)                                             \
  PRECICE_CHECK(!_meshLock.check(name),                                               \
                "This participant attempted to modify the Mesh \"{}\" while locked. " \
                "Mesh modification is only allowed before calling initialize().",     \
                name);

/** Validates a given meshID
 * This macros creates the "id" in a local scope and provides it to the called implementation.
 */
#define PRECICE_VALIDATE_MESH_NAME(name)  \
  do {                                    \
    PRECICE_VALIDATE_MESH_NAME_IMPL(name) \
  } while (false)

/** Validates a given meshID and checks if the mesh is used by the current participant
 * This macros creates the "id" in a local scope and provides it to the called implementation.
 */
#define PRECICE_REQUIRE_MESH_USE(name)  \
  do {                                  \
    PRECICE_REQUIRE_MESH_USE_IMPL(name) \
  } while (false)

/** Validates a given meshID and checks if the mesh is provided by the current participant
 * This macros creates the "id" in a local scope and provides it to the called implementation.
 */
#define PRECICE_REQUIRE_MESH_PROVIDE(name)  \
  do {                                      \
    PRECICE_REQUIRE_MESH_PROVIDE_IMPL(name) \
  } while (false)

/** Validates a given meshID, checks if the mesh is provided by the current participant and unlocked
 * This macros creates the "id" in a local scope and provides it to the called implementation.
 */
#define PRECICE_REQUIRE_MESH_MODIFY(name)  \
  do {                                     \
    PRECICE_REQUIRE_MESH_MODIFY_IMPL(name) \
  } while (false)

//
// DATA VALIDATION
//

/** Implementation of PRECICE_VALIDATE_DATA_ID()
 *
 * @attention Do not use this macro directly!
 */
#define PRECICE_VALIDATE_DATA_NAME_IMPL(mesh, data) \
  PRECICE_CHECK(_accessor->hasData(mesh, data),     \
                "The given data \"{}\" is unknown to preCICE mesh \"{}\".{}", data, mesh, _accessor->hintForMeshData(mesh, data));

/** Implementation of PRECICE_REQUIRE_DATA_READ()
 *
 * @attention Do not use this macro directly!
 */
#define PRECICE_REQUIRE_DATA_READ_IMPL(mesh, data)                                                                             \
  PRECICE_VALIDATE_DATA_NAME_IMPL(mesh, data)                                                                                  \
  PRECICE_CHECK((_accessor->isDataRead(mesh, data)),                                                                           \
                "This participant does not use Data \"{0}\", but attempted to read it. "                                       \
                "Please extend the configuration of participant \"{1}\" by defining <read-data mesh=\"{2}\" name=\"{0}\" />.", \
                data, _accessorName, mesh);

/** Implementation of PRECICE_REQUIRE_DATA_WRITE()
 *
 * @attention Do not use this macro directly!
 */
#define PRECICE_REQUIRE_DATA_WRITE_IMPL(mesh, data)                                                                             \
  PRECICE_VALIDATE_DATA_NAME_IMPL(mesh, data)                                                                                   \
  PRECICE_CHECK((_accessor->isDataWrite(mesh, data)),                                                                           \
                "This participant does not use Data \"{0}\", but attempted to write it. "                                       \
                "Please extend the configuration of participant \"{1}\" by defining <write-data mesh=\"{2}\" name=\"{0}\" />.", \
                data, _accessorName, mesh);

/** Validates a given dataID
 * This macros creates the "id" in a local scope and provides it to the called implementation.
 */
#define PRECICE_VALIDATE_DATA_NAME(mesh, data)  \
  do {                                          \
    PRECICE_VALIDATE_DATA_NAME_IMPL(mesh, data) \
  } while (false)

/** Validates a given dataID and checks for read access
 * This macros creates the "id" in a local scope and provides it to the called implementation.
 */
#define PRECICE_REQUIRE_DATA_READ(mesh, data)  \
  do {                                         \
    PRECICE_REQUIRE_DATA_READ_IMPL(mesh, data) \
  } while (false)

/** Validates a given dataID and checks for write access
 * This macros creates the "id" in a local scope and provides it to the called implementation.
 */
#define PRECICE_REQUIRE_DATA_WRITE(mesh, data)  \
  do {                                          \
    PRECICE_REQUIRE_DATA_WRITE_IMPL(mesh, data) \
  } while (false)

//
// DATA VALUE VALIDATION
//
#ifdef NDEBUG

#define PRECICE_VALIDATE_DATA(data, size) \
  {                                       \
  }
#else // NDEBUG

#define PRECICE_VALIDATE_DATA(data, size) \
  PRECICE_CHECK(std::all_of(data, data + size, [](double val) { return std::isfinite(val); }), "One of the given data values is either plus or minus infinity or NaN.");

#endif

#define PRECICE_EXPERIMENTAL_API()                                                                                                                    \
  PRECICE_CHECK(_allowsExperimental, "You called the API function \"{}\", which is part of the experimental API. "                                    \
                                     "You may unlock the full API by specifying <solver-interface experimental=\"true\" ... > in the configuration. " \
                                     "Please be aware that experimental features may change in any future version (even minor or bugfix).",           \
                __func__)
