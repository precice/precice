#pragma once

/**@file
 * This file contains verification macros to harden the SolverInterface against misuse and misconfiguration.
 *
 * \attention Only include this file from SolverInterfaceImpl.cpp
 */

//
// MESH VALIDATION
//

/** Implementation of PRECICE_VALIDATE_MESH_ID()
 *
 * @attention Do not use this macro directly!
 */
#include "utils/stacktrace.hpp"
#define PRECICE_VALIDATE_MESH_ID_IMPL(id)                          \
  PRECICE_CHECK(_accessor->hasMesh(id),                            \
                "The given Mesh ID \"{}\" is unknown to preCICE.", \
                id);

/** Implementation of PRECICE_REQUIRE_MESH_USE()
 *
 * @attention Do not use this macro directly!
 */
#define PRECICE_REQUIRE_MESH_USE_IMPL(id)                                                              \
  PRECICE_VALIDATE_MESH_ID_IMPL(id)                                                                    \
  PRECICE_CHECK(_accessor->isMeshUsed(id),                                                             \
                "This participant does not use the mesh \"{0}\", but attempted to access it. "         \
                "Please define <use-mesh name=\"{0}\" /> in the configuration of participant \" {1}.", \
                _accessor->getMeshName(id), _accessorName);

/** Implementation of PRECICE_REQUIRE_MESH_PROVIDE()
 *
 * @attention Do not use this macro directly!
 */
#define PRECICE_REQUIRE_MESH_PROVIDE_IMPL(id)                                                           \
  PRECICE_REQUIRE_MESH_USE_IMPL(id)                                                                     \
  PRECICE_CHECK(_accessor->isMeshProvided(id),                                                          \
                "This participant does not provide Mesh \"{0}\", but attempted to modify it. "          \
                "Please extend the use-mesh tag as follows <use-mesh name=\"{0}\" provide=\"yes\" />.", \
                _accessor->getMeshName(id));

/** Implementation of PRECICE_REQUIRE_MESH_MODIFY()
 *
 * @attention Do not use this macro directly!
 */
#define PRECICE_REQUIRE_MESH_MODIFY_IMPL(id)                                          \
  PRECICE_REQUIRE_MESH_PROVIDE_IMPL(id)                                               \
  PRECICE_CHECK(!_meshLock.check(id),                                                 \
                "This participant attempted to modify the Mesh \"{}\" while locked. " \
                "Mesh modification is only allowed before calling initialize().",     \
                _accessor->getMeshName(id));

/** Validates a given meshID
 * This macros creates the "id" in a local scope and provides it to the called implementation.
 */
#define PRECICE_VALIDATE_MESH_ID(meshID) \
  do {                                   \
    const auto id = (meshID);            \
    PRECICE_VALIDATE_MESH_ID_IMPL(id)    \
  } while (false)

/** Validates a given meshID and checks if the mesh is used by the current participant
 * This macros creates the "id" in a local scope and provides it to the called implementation.
 */
#define PRECICE_REQUIRE_MESH_USE(meshID) \
  do {                                   \
    const auto id = (meshID);            \
    PRECICE_REQUIRE_MESH_USE_IMPL(id)    \
  } while (false)

/** Validates a given meshID and checks if the mesh is provided by the current participant
 * This macros creates the "id" in a local scope and provides it to the called implementation.
 */
#define PRECICE_REQUIRE_MESH_PROVIDE(meshID) \
  do {                                       \
    const auto id = (meshID);                \
    PRECICE_REQUIRE_MESH_PROVIDE_IMPL(id)    \
  } while (false)

/** Validates a given meshID, checks if the mesh is provided by the current participant and unlocked
 * This macros creates the "id" in a local scope and provides it to the called implementation.
 */
#define PRECICE_REQUIRE_MESH_MODIFY(meshID) \
  do {                                      \
    const auto id = (meshID);               \
    PRECICE_REQUIRE_MESH_MODIFY_IMPL(id)    \
  } while (false)

//
// DATA VALIDATION
//

/** Implementation of PRECICE_VALIDATE_DATA_ID()
 *
 * @attention Do not use this macro directly!
 */
#define PRECICE_VALIDATE_DATA_ID_IMPL(id) \
  PRECICE_CHECK(_accessor->hasData(id),   \
                "The given Data ID \"{}\" is unknown to preCICE.", id);

/** Implementation of PRECICE_REQUIRE_DATA_READ()
 *
 * @attention Do not use this macro directly!
 */
#define PRECICE_REQUIRE_DATA_READ_IMPL(id)                                                                                      \
  PRECICE_VALIDATE_DATA_ID_IMPL(id)                                                                                             \
  PRECICE_CHECK((_accessor->isDataUsed(id) && _accessor->isDataRead(id)),                                                       \
                "This participant does not use Data \"{0}\", but attempted to read it. "                                        \
                "Please extend the configuarion of partiticipant \"{1}\" by defining <read-data mesh=\"{0}\" name=\"{2}\" />.", \
                _accessor->getDataName(id), _accessorName, _accessor->getMeshNameFromData(id));

/** Implementation of PRECICE_REQUIRE_DATA_WRITE()
 *
 * @attention Do not use this macro directly!
 */
#define PRECICE_REQUIRE_DATA_WRITE_IMPL(id)                                                                                      \
  PRECICE_VALIDATE_DATA_ID_IMPL(id)                                                                                              \
  PRECICE_CHECK((_accessor->isDataUsed(id) && _accessor->isDataWrite(id)),                                                       \
                "This participant does not use Data \"{0}\", but attempted to write it. "                                        \
                "Please extend the configuarion of partiticipant \"{1}\" by defining <write-data mesh=\"{0}\" name=\"{2}\" />.", \
                _accessor->getDataName(id), _accessorName, _accessor->getMeshNameFromData(id));

/** Validates a given dataID
 * This macros creates the "id" in a local scope and provides it to the called implementation.
 */
#define PRECICE_VALIDATE_DATA_ID(dataID) \
  do {                                   \
    const auto id = (dataID);            \
    PRECICE_VALIDATE_DATA_ID_IMPL(id)    \
  } while (false)

/** Validates a given dataID and checks for read access
 * This macros creates the "id" in a local scope and provides it to the called implementation.
 */
#define PRECICE_REQUIRE_DATA_READ(dataID) \
  do {                                    \
    const auto id = (dataID);             \
    PRECICE_REQUIRE_DATA_READ_IMPL(id)    \
  } while (false)

/** Validates a given dataID and checks for write access
 * This macros creates the "id" in a local scope and provides it to the called implementation.
 */
#define PRECICE_REQUIRE_DATA_WRITE(dataID) \
  do {                                     \
    const auto id = (dataID);              \
    PRECICE_REQUIRE_DATA_WRITE_IMPL(id)    \
  } while (false)

//
// DATA VALUE VALIDATION
//
#ifdef NDEBUG

#define PRECICE_VALIDATE_DATA(data, size) \
  {                                       \
  }
#else //NDEBUG

#define PRECICE_VALIDATE_DATA(data, size) \
  PRECICE_CHECK(std::all_of(data, data + size, [](double val) { return std::isfinite(val); }), "One of the given data values is either plus or minus infinity or NaN.");

#endif
