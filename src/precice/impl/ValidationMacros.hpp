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
#define PRECICE_VALIDATE_MESH_ID_IMPL(id) \
  PRECICE_CHECK(_dataIDs.find(id) != _dataIDs.end(), "The given Mesh ID \"" << id << "\" is unknown to preCICE.");

/** Implementation of PRECICE_REQUIRE_MESH_USE()
 *
 * @attention Do not use this macro directly!
 */
#define PRECICE_REQUIRE_MESH_USE_IMPL(id)            \
  PRECICE_VALIDATE_MESH_ID_IMPL(id)                  \
  MeshContext &context = _accessor->meshContext(id); \
  PRECICE_CHECK(_accessor->isMeshUsed(id), "This participant does not use the mesh \"" << context.mesh->getName() << "\", but attempted to access it. Please define <use-mesh name=\"" << context.mesh->getName() << "\" /> in the configuration of participant \"" << _accessorName << '"');

/** Implementation of PRECICE_REQUIRE_MESH_PROVIDE()
 *
 * @attention Do not use this macro directly!
 */
#define PRECICE_REQUIRE_MESH_PROVIDE_IMPL(id) \
  PRECICE_REQUIRE_MESH_USE_IMPL(id)           \
  PRECICE_CHECK(context.provideMesh, "This participant does not provide Mesh \"" << context.mesh->getName() << "\", but attempted to modify it. Please extend the use-mesh tag as follows <use-mesh name=\"" << context.mesh->getName() << "\" provide=\"yes\" />.");

/** Implementation of PRECICE_REQUIRE_MESH_MODIFY()
 *
 * @attention Do not use this macro directly!
 */
#define PRECICE_REQUIRE_MESH_MODIFY_IMPL(id) \
  PRECICE_REQUIRE_MESH_PROVIDE_IMPL(id)      \
  PRECICE_CHECK(!_meshLock.check(meshID), "This participant attempted to modify the Mesh \"" << context.mesh->getName() << "\" while locked. Mesh modification is only allowed before calling initialize().");

/// Validates a given meshID
#define PRECICE_VALIDATE_MESH_ID(meshID) \
  do {                                   \
    const auto id = (meshID);            \
    PRECICE_VALIDATE_MESH_ID_IMPL(id)    \
  } while (false)

/// Validates a given meshID and checks if the mesh is used by the current participant
#define PRECICE_REQUIRE_MESH_USE(meshID) \
  do {                                   \
    const auto id = (meshID);            \
    PRECICE_REQUIRE_MESH_USE_IMPL(id)    \
  } while (false)

/// Validates a given meshID and checks if the mesh is provided by the current participant
#define PRECICE_REQUIRE_MESH_PROVIDE(meshID) \
  do {                                       \
    const auto id = (meshID);                \
    PRECICE_REQUIRE_MESH_PROVIDE_IMPL(id)    \
  } while (false)

/// Validates a given meshID and checks if the mesh is provided by the current participant and unlocked
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
#define PRECICE_VALIDATE_DATA_ID_IMPL(id)                                                                                                           \
  PRECICE_CHECK(std::any_of(_dataIDs.begin(), _dataIDs.end(), [id](const typename decltype(_dataIDs)::value_type &meshkv) {                         \
                  return std::any_of(meshkv.second.begin(), meshkv.second.end(), [id](const typename decltype(meshkv.second)::value_type &datakv) { \
                    return datakv.second == id;                                                                                                     \
                  });                                                                                                                               \
                }),                                                                                                                                 \
                "The given Data ID \"" << id << "\" is unknown to preCICE.");

/** Implementation of PRECICE_REQUIRE_DATA_READ()
 *
 * @attention Do not use this macro directly!
 */
#define PRECICE_REQUIRE_DATA_READ_IMPL(id)                                                                                                                        \
  PRECICE_VALIDATE_DATA_ID_IMPL(id)                                                                                                                               \
  DataContext &context = _accessor->dataContext(id);                                                                                                              \
  PRECICE_CHECK((_accessor->isDataUsed(id) && _accessor->isDataRead(id)),                                                                                         \
                "This participant does not use Data \"" << context.getName() << "\", but attempted to read it. Please extend the configuarion of partiticpant \"" \
                                                        << _accessorName << "\" by defining <read-data mesh=\"" << context.mesh->getName() << "\" name=\"" << context.getName() << "\" />.");

/** Implementation of PRECICE_REQUIRE_DATA_WRITE()
 *
 * @attention Do not use this macro directly!
 */
#define PRECICE_REQUIRE_DATA_WRITE_IMPL(id)                                                                                                                        \
  PRECICE_VALIDATE_DATA_ID_IMPL(id)                                                                                                                                \
  DataContext &context = _accessor->dataContext(id);                                                                                                               \
  PRECICE_CHECK((_accessor->isDataUsed(id) && _accessor->isDataWrite(id)),                                                                                         \
                "This participant does not use Data \"" << context.getName() << "\", but attempted to write it. Please extend the configuarion of partiticpant \"" \
                                                        << _accessorName << "\" by defining <write-data mesh=\"" << context.mesh->getName() << "\" name=\"" << context.getName() << "\" />.");

/// Validates a given dataID
#define PRECICE_VALIDATE_DATA_ID(dataID) \
  do {                                   \
    const auto id = (dataID);            \
    PRECICE_VALIDATE_DATA_ID_IMPL(id)    \
  } while (false)

/// Validates a dataID and checks for read access
#define PRECICE_REQUIRE_DATA_READ(dataID) \
  do {                                    \
    const auto id = (dataID);             \
    PRECICE_REQUIRE_DATA_READ_IMPL(id)    \
  } while (false)

/// Validates a dataID and checks for write access
#define PRECICE_REQUIRE_DATA_WRITE(dataID) \
  do {                                     \
    const auto id = (dataID);              \
    PRECICE_REQUIRE_DATA_WRITE_IMPL(id)    \
  } while (false)
