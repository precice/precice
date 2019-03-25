#pragma once

/**@file
 * This file contains verification macros to harden the SolverInterface against misuse and misconfiguration.
 *
 * \attention Only include this file from SolverInterfaceImpl.cpp
 */

//
// MESH VALIDATION 
//

/// Implementation of PRECICE_VALIDATE_MESH_ID
#define PRECICE_VALIDATE_MESH_ID_IMPL(id) \
    CHECK(_dataIDs.find(id) != _dataIDs.end(), "There is no Mesh with ID:" << id);

/// Implementation of PRECICE_REQUIRE_MESH_USE
#define PRECICE_REQUIRE_MESH_USE_IMPL(id) \
    PRECICE_VALIDATE_MESH_ID_IMPL(id) \
    CHECK(_accessor->isMeshUsed(id), "Mesh is not used by this participant! ID:" << id);

/// Implementation of PRECICE_REQUIRE_MESH_PROVIDE
#define PRECICE_REQUIRE_MESH_PROVIDE_IMPL(id) \
    PRECICE_REQUIRE_MESH_USE_IMPL(id) \
    MeshContext& context = _accessor->meshContext(id); \
    CHECK(context.provideMesh, "Can only write to provided meshes! " << id);

/// Validates a given meshID
#define PRECICE_VALIDATE_MESH_ID(meshID) \
    do { \
        const auto id = (meshID); \
        PRECICE_VALIDATE_MESH_ID_IMPL(id) \
    } while(false)

/// Validates a given meshID and checks if the mesh is used by the current participant
#define PRECICE_REQUIRE_MESH_USE(meshID) \
    do { \
        const auto id = (meshID); \
        PRECICE_REQUIRE_MESH_USE_IMPL(id) \
    } while(false)

/// Validates a given meshID and checks if the mesh is provided by the current participant
#define PRECICE_REQUIRE_MESH_PROVIDE(meshID) \
    do { \
        const auto id = (meshID); \
        PRECICE_REQUIRE_MESH_PROVIDE_IMPL(id) \
    } while(false)

//
// DATA VALIDATION 
//

/// Implementation of PRECICE_VALIDATE_DATA_ID
#define PRECICE_VALIDATE_DATA_ID_IMPL(id) \
    CHECK(std::any_of(_dataIDs.begin(), _dataIDs.end(), [id](const typename decltype(_dataIDs)::value_type & meshkv){ \
                return std::any_of(meshkv.second.begin(), meshkv.second.end(), [id](const typename decltype(meshkv.second)::value_type& datakv){ \
                        return datakv.second == id; \
                        }); \
                }), \
            "There is no Data with ID " << id); \

/// Implementation of PRECICE_REQUIRE_DATA_READ
#define PRECICE_REQUIRE_DATA_READ_IMPL(id) \
    PRECICE_VALIDATE_DATA_ID_IMPL(id) \
    CHECK(_accessor->isDataUsed(id), "Data is not used by this participant! " << id); \
    CHECK(_accessor->isDataRead(id), "Data is not marked as read! " << id); \

/// Implementation of PRECICE_REQUIRE_DATA_WRITE
#define PRECICE_REQUIRE_DATA_WRITE_IMPL(id) \
    PRECICE_VALIDATE_DATA_ID_IMPL(id) \
    CHECK(_accessor->isDataUsed(id), "Data is not used by this participant! " << id); \
    CHECK(_accessor->isDataWrite(id), "Data is not marked as write! " << id); \

/// Validates a given dataID
#define PRECICE_VALIDATE_DATA_ID(dataID) \
    do { \
        const auto id = (dataID); \
        PRECICE_VALIDATE_DATA_ID_IMPL(id) \
    } while(false)

/// Validates a dataID and checks for read access
#define PRECICE_REQUIRE_DATA_READ(dataID) \
    do { \
        const auto id = (dataID); \
        PRECICE_REQUIRE_DATA_READ_IMPL(id) \
    } while(false)

/// Validates a dataID and checks for write access
#define PRECICE_REQUIRE_DATA_WRITE(dataID) \
    do { \
        const auto id = (dataID); \
        PRECICE_REQUIRE_DATA_WRITE_IMPL(id) \
    } while(false)
