#pragma once

/** @file
 * This file contains a Fortran 77 compatible interface written in C/C++.
 *
 *
 * Every method has a Fortran syntax equivalent in the method comment, and a
 * listing for input and output variables. A variable can be input and output
 * at the same time.
 */

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Fortran syntax:
 * precicef_create(
 *   CHARACTER participantName(*),
 *   CHARACTER configFileName(*),
 *   INTEGER   solverProcessIndex,
 *   INTEGER   solverProcessSize )
 *
 * IN:  participantName, configFileName, solverProcessIndex, solverProcessSize
 * OUT: -
 *
 * @copydoc precice::SolverInterface::SolverInterface()
 *
 */
void precicef_create_(
    const char *participantName,
    const char *configFileName,
    const int * solverProcessIndex,
    const int * solverProcessSize,
    int         lengthAccessorName,
    int         lengthConfigFileName);

/**
 * Fortran syntax:
 * precicef_initialize( DOUBLE PRECISION timstepLengthLimit )
 *
 * IN:  -
 * OUT: timestepLengthLimit
 *
 * @copydoc precice::SolverInterface::initialize()
 *
 */
void precicef_initialize_(double *timestepLengthLimit);

/**
 * Fortran syntax:
 * precicef_intialize_data()
 *
 * IN: -
 * OUT: -
 *
 * @copydoc precice::SolverInterface::initializeData()
 *
 */
void precicef_initialize_data_();

/**
 * Fortran syntax:
 * precicef_advance( DOUBLE PRECISION timstepLengthLimit )
 *
 * IN:  timestepLengthLimit
 * OUT: timestepLengthLimit
 *
 * @copydoc precice::SolverInterface::advance()
 *
 */
void precicef_advance_(double *timestepLengthLimit);

/**
 * Fortran syntax:
 * precicef_finalize();
 *
 * @copydoc precice::SolverInterface::finalize()
 *
 */
void precicef_finalize_();

/**
 * Fortran syntax:
 * precicef_get_dims( INTEGER dimensions )
 *
 * IN:  -
 * OUT: dimensions
 *
 * @copydoc precice::SolverInterface::getDimensions()
 *
 */
void precicef_get_dims_(int *dimensions);

/**
 * @deprecated Forwards to precicef_is_coupling_ongoing_
 *
 * Fortran syntax:
 * precicef_ongoing( INTEGER isOngoing )
 *
 * IN:  -
 * OUT: isOngoing(1:true, 0:false)
 *
 * @copydoc precice::SolverInterface::isOngoing()
 *
 */
[[deprecated("Use precicef_is_coupling_ongoing_() instead.")]] void precicef_ongoing_(int *isOngoing);

/**
 * Fortran syntax:
 * precicef_is_coupling_ongoing( INTEGER isOngoing )
 *
 * IN:  -
 * OUT: isOngoing(1:true, 0:false)
 *
 * @copydoc precice::SolverInterface::isCouplingOngoing()
 *
 */
void precicef_is_coupling_ongoing_(int *isOngoing);

/**
 * @deprecated Forwards to precicef_is_write_data_required_
 *
 * Fortran syntax:
 * precicef_write_data_required(
 *  DOUBLE PRECISION computedTimestepLength,
 *  INTEGER          isRequired )
 *
 * IN:  computedTimestepLength
 * OUT: isRequired(1:true, 0:false)
 *
 * @copydoc precice::SolverInterface::isWriteDataRequired()
 *
 */
[[deprecated("Use precicef_is_write_data_required_(...) with the same arguments instead.")]] void precicef_write_data_required_(
    const double *computedTimestepLength,
    int *         isRequired);

/**
 * Fortran syntax:
 * precicef_is_write_data_required(
 *  DOUBLE PRECISION computedTimestepLength,
 *  INTEGER          isRequired )
 *
 * IN:  computedTimestepLength
 * OUT: isRequired(1:true, 0:false)
 *
 * @copydoc precice::SolverInterface::isWriteDataRequired()
 *
 */
void precicef_is_write_data_required_(
    const double *computedTimestepLength,
    int *         isRequired);

/**
 * @deprecated Forwards to precicef_is_read_data_available_
 *
 * Fortran syntax:
 * precicef_read_data_available( INTEGER isAvailable );
 *
 * IN:  -
 * OUT: isAvailable(1:true, 0:false)
 *
 * @copydoc precice::SolverInterface::isReadDataAvailable()
 *
 */
[[deprecated("Use precicef_is_read_data_available_() instead.")]] void precicef_read_data_available_(int *isAvailable);

/**
 * Fortran syntax:
 * precicef_is_read_data_available( INTEGER isAvailable );
 *
 * IN:  -
 * OUT: isAvailable(1:true, 0:false)
 *
 * @copydoc precice::SolverInterface::isReadDataAvailable()
 *
 */
void precicef_is_read_data_available_(int *isAvailable);

/**
 * Fortran syntax:
 * precicef_is_time_window_complete( INTEGER isComplete );
 *
 * IN:  -
 * OUT: isComplete(1:true, 0:false)
 *
 * @copydoc precice::SolverInterface::isTimeWindowComplete()
 *
 */
void precicef_is_time_window_complete_(int *isComplete);

/**
 * Fortran syntax:
 * precicef_has_to_evaluate_surrogate_model( INTEGER hasToEvaluate );
 *
 * IN:  -
 * OUT: hasToEvaluate(1:true, 0:false)
 *
 * @copydoc precice::SolverInterface::hasToEvaluateSurrogateModel()
 *
 */
void precicef_has_to_evaluate_surrogate_model_(int *hasToEvaluate);

/**
 * Fortran syntax:
 * precicef_has_to_evaluate_fine_model( INTEGER hasToEvaluate );
 *
 * IN:  -
 * OUT: hasToEvaluate(1:true, 0:false)
 *
 * @copydoc precice::SolverInterface::hasToEvaluateFineModel()
 *
 */
void precicef_has_to_evaluate_fine_model_(int *hasToEvaluate);

/**
 * @deprecated Forwards to precicef_is_action_required_
 *
 * Fortran syntax:
 * precicef_action_required(
 *   CHARACTER action(*),
 *   INTEGER   isRequired )
 *
 * IN:  action
 * OUT: isRequired(1:true, 0:false)
 *
 * @copydoc precice::SolverInterface::isActionRequired()
 *
 */
[[deprecated("Use precicef_is_action_required_(...) with the same arguments instead.")]] void precicef_action_required_(
    const char *action,
    int *       isRequired,
    int         lengthAction);

/**
 * Fortran syntax:
 * precicef_is_action_required(
 *   CHARACTER action(*),
 *   INTEGER   isRequired )
 *
 * IN:  action
 * OUT: isRequired(1:true, 0:false)
 *
 * @copydoc precice::SolverInterface::isActionRequired()
 *
 */
void precicef_is_action_required_(
    const char *action,
    int *       isRequired,
    int         lengthAction);

/**
 * Fortran syntax:
 * precicef_mark_action_fulfilled( CHARACTER action(*) )
 *
 * IN:  action
 * OUT: -
 *
 * @copydoc precice::SolverInterface::markActionFulfilled()
 *
 */
void precicef_mark_action_fulfilled_(
    const char *action,
    int         lengthAction);

/**
 * Fortran syntax:
 * precicef_has_mesh(
 *   CHARACTER meshName(*),
 *   INTEGER   hasMesh )
 *
 * IN:  meshName
 * OUT: hasMesh(1:true, 0:false)
 *
 * @copydoc precice::SolverInterface::hasMesh()
 *
 */
void precicef_has_mesh_(
    const char *meshName,
    int *       hasMesh,
    int         lengthMeshName);

/**
 * Fortran syntax:
 * precicef_get_mesh_id(
 *   CHARACTER meshName(*),
 *   INTEGER   meshID )
 *
 * IN:  meshName
 * OUT: meshID
 *
 * @copydoc precice::SolverInterface::getMeshID()
 *
 */
void precicef_get_mesh_id_(
    const char *meshName,
    int *       meshID,
    int         lengthMeshName);

/**
 * Fortran syntax:
 * precicef_has_data(
 *   CHARACTER dataName(*),
 *   INTEGER   meshID,
 *   INTEGER   hasData)
 *
 * IN:  dataName
 * IN:  meshID
 * OUT: hasData(1:true, 0:false)
 *
 * @copydoc precice::SolverInterface::hasData()
 *
 */
void precicef_has_data_(
    const char *dataName,
    const int * meshID,
    int *       hasData,
    int         lengthDataName);

/**
 * The given name (dataName) has to be one of the names specified in the
 * configuration file. The data id obtained can be used to read and write
 * data to and from the coupling mesh.
 *
 * Fortran syntax:
 * precicef_get_data_id(
 *   CHARACTER dataName(*),
 *   INTEGER   meshID,
 *   INTEGER   dataID,
 *   INTEGER   lengthDataName)
 *
 * IN:  dataName
 * IN:  meshID
 * OUT: dataID
 *
 * @copydoc precice::SolverInterface::getDataID()
 *
 */
void precicef_get_data_id_(
    const char *dataName,
    const int * meshID,
    int *       dataID,
    int         lengthDataName);

/**
 * Fortran syntax:
 * precicef_has_data(
 *   INTEGER   meshID
 *   INTEGER   required)
 *
 * IN:  meshID
 * OUT: required(1:true, 0:false)
 *
 * @copydoc precice::SolverInterface::isMeshConnectivityRequired()
 */
void precicef_is_mesh_connectivity_required_(
    const int *meshID,
    int *      required);

/**
 * Fortran syntax:
 * precicef_set_vertex(
 *   INTEGER          meshID,
 *   DOUBLE PRECISION position(dim),
 *   INTEGER          vertexID )
 *
 * IN:  meshID, position
 * OUT: vertexID
 *
 * @copydoc precice::SolverInterface::setMeshVertex()
 *
 */
void precicef_set_vertex_(
    const int *   meshID,
    const double *position,
    int *         vertexID);

/**
 * Fortran syntax:
 * precicef_get_mesh_vertex_size(
 *   INTEGER meshID,
 *   INTEGER meshSize )
 *
 * IN:  meshID
 * OUT: meshSize
 *
 * @copydoc precice::SolverInterface::getMeshVertexSize()
 *
 */
void precicef_get_mesh_vertex_size_(
    const int *meshID,
    int *      meshSize);

/**
 * Fortran syntax:
 * precicef_set_vertices(
 *   INTEGER          meshID,
 *   INTEGER          size,
 *   DOUBLE PRECISION positions(dim*size),
 *   INTEGER          positionIDs(size) )
 *
 * IN:  meshID, size, positions
 * OUT: positionIDs
 *
 * @copydoc precice::SolverInterface::setMeshVertices()
 *
 */
void precicef_set_vertices_(
    const int *meshID,
    const int *size,
    double *   positions,
    int *      positionIDs);

/**
 * Fortran syntax:
 * precicef_get_vertices(
 *   INTEGER          meshID,
 *   INTEGER          size,
 *   INTEGER          ids(size)
 *   DOUBLE PRECISION positions(dim*size))
 *
 * IN:  meshID, size, ids
 * OUT: positions
 *
 * @copydoc precice::SolverInterface::getMeshVertices()
 *
 */
void precicef_get_vertices_(
    const int *meshID,
    const int *size,
    int *      ids,
    double *   positions);

/**
 * Fortran syntax:
 * precicef_get_vertices(
 *   INTEGER          meshID,
 *   INTEGER          size,
 *   DOUBLE PRECISION positions(dim*size),
 *   INTEGER          ids(size))
 *
 * IN:  meshID, size, positions
 * OUT: ids
 *
 * @copydoc precice::SolverInterface::getMeshVertexIDsFromPositions()
 *
 */
void precicef_get_vertex_ids_from_positions_(
    const int *meshID,
    const int *size,
    double *   positions,
    int *      ids);

/**
 * Fortran syntax:
 * precicef_set_edge(
 *   INTEGER meshID,
 *   INTEGER firstVertexID,
 *   INTEGER secondVertexID,
 *   INTEGER edgeID )
 *
 * IN:  meshID, firstVertexID, secondVertexID
 * OUT: edgeID
 *
 * @copydoc precice::SolverInterface::setMeshEdge()
 *
 */
void precicef_set_edge_(
    const int *meshID,
    const int *firstVertexID,
    const int *secondVertexID,
    int *      edgeID);

/**
 * Fortran syntax:
 * precicef_set_triangle(
 *   INTEGER meshID,
 *   INTEGER firstEdgeID,
 *   INTEGER secondEdgeID,
 *   INTEGER thirdEdgeID )
 *
 * IN:  meshID, firstEdgeID, secondEdgeID, thirdEdgeID
 * OUT: -
 *
 * @copydoc precice::SolverInterface::setMeshTriangle()
 *
 */
void precicef_set_triangle_(
    const int *meshID,
    const int *firstEdgeID,
    const int *secondEdgeID,
    const int *thirdEdgeID);

/**
 * Fortran syntax:
 * precicef_set_triangle_we(
 *   INTEGER meshID,
 *   INTEGER firstVertexID,
 *   INTEGER secondVertexID,
 *   INTEGER thirdVertexID )
 *
 * IN:  meshID, firstVertexID, secondVertexID, thirdVertexID
 * OUT: -
 *
 * @copydoc precice::SolverInterface::setMeshTriangleWithEdges()
 *
 */
void precicef_set_triangle_we_(
    const int *meshID,
    const int *firstVertexID,
    const int *secondVertexID,
    const int *thirdVertexID);

/**
 * Fortran syntax:
 * precicef_set_quad(
 *   INTEGER meshID,
 *   INTEGER firstEdgeID,
 *   INTEGER secondEdgeID,
 *   INTEGER thirdEdgeID,
 *   INTEGER fourthEdgeID )
 *
 * IN:  meshID, firstEdgeID, secondEdgeID, thirdEdgeID, fourthEdgeID
 * OUT: -
 *
 * @copydoc precice::SolverInterface::setMeshQuad()
 *
 */
void precicef_set_quad_(
    const int *meshID,
    const int *firstEdgeID,
    const int *secondEdgeID,
    const int *thirdEdgeID,
    const int *fourthEdgeID);

/**
 * Fortran syntax:
 * precicef_set_quad_we(
 *   INTEGER meshID,
 *   INTEGER firstVertexID,
 *   INTEGER secondVertexID,
 *   INTEGER thirdVertexID,
 *   INTEGER fourthVertexID )
 *
 * IN:  meshID, firstVertexID, secondVertexID, thirdVertexID, fourthVertexID
 * OUT: -
 *
 * @copydoc precice::SolverInterface::setMeshQuadWithEdges()
 *
 */
void precicef_set_quad_we_(
    const int *meshID,
    const int *firstVertexID,
    const int *secondVertexID,
    const int *thirdVertexID,
    const int *fourthVertexID);

/**
 * Fortran syntax:
 * precicef_write_bvdata(
 *   INTEGER dataID,
 *   INTEGER size,
 *   INTEGER valueIndices,
 *   DOUBLE PRECISION values(dim*size) )
 *
 * IN:  dataID, size, valueIndices, values
 * OUT: -
 *
 * @copydoc precice::SolverInterface::writeBlockVectorData
 *
 */
void precicef_write_bvdata_(
    const int *dataID,
    const int *size,
    int *      valueIndices,
    double *   values);

/**
 * Fortran syntax:
 * precicef_write_vdata(
 *   INTEGER dataID,
 *   INTEGER valueIndex,
 *   DOUBLE PRECISION dataValue(dim) )
 *
 * IN:  dataID, valueIndex, dataValue
 * OUT: -
 *
 * @copydoc precice::SolverInterface::writeVectorData
 *
 */
void precicef_write_vdata_(
    const int *   dataID,
    const int *   valueIndex,
    const double *dataValue);

/**
 * Fortran syntax:
 * precicef_write_bsdata(
 *   INTEGER dataID,
 *   INTEGER size,
 *   INTEGER valueIndices,
 *   DOUBLE PRECISION values(size) )
 *
 * IN:  dataID, size, valueIndices, values
 * OUT: -
 *
 * @copydoc precice::SolverInterface::writeBlockScalarData
 *
 */
void precicef_write_bsdata_(
    const int *dataID,
    const int *size,
    int *      valueIndices,
    double *   values);

/**
 * Fortran syntax:
 * precicef_write_sdata(
 *   INTEGER dataID,
 *   INTEGER valueIndex,
 *   DOUBLE PRECISION dataValue )
 *
 * IN:  dataID, valueIndex, dataValue
 * OUT: -
 *
 * @copydoc precice::SolverInterface::writeScalarData
 *
 */
void precicef_write_sdata_(
    const int *   dataID,
    const int *   valueIndex,
    const double *dataValue);

/**
 * Fortran syntax:
 * precicef_read_bvdata(
 *   INTEGER dataID,
 *   INTEGER size,
 *   INTEGER valueIndices,
 *   DOUBLE PRECISION values(dim*size) )
 *
 * IN:  dataID, size, valueIndices
 * OUT: values
 *
 * @copydoc precice::SolverInterface::readBlockVectorData
 *
 */
void precicef_read_bvdata_(
    const int *dataID,
    const int *size,
    int *      valueIndices,
    double *   values);

/**
 * Fortran syntax:
 * precicef_read_vdata(
 *   INTEGER dataID,
 *   INTEGER valueIndex,
 *   DOUBLE PRECISION dataValue(dim) )
 *
 * IN:  dataID, valueIndex
 * OUT: dataValue
 *
 * @copydoc precice::SolverInterface::readVectorData
 *
 */
void precicef_read_vdata_(
    const int *dataID,
    const int *valueIndex,
    double *   dataValue);

/**
 * Fortran syntax:
 * precicef_read_bsdata(
 *   INTEGER dataID,
 *   INTEGER size,
 *   INTEGER valueIndices,
 *   DOUBLE PRECISION values(size) )
 *
 * IN:  dataID, size, valueIndices
 * OUT: values
 *
 * @copydoc precice::SolverInterface::readBlockScalarData
 *
 */
void precicef_read_bsdata_(
    const int *dataID,
    const int *size,
    int *      valueIndices,
    double *   values);

/**
 * Fortran syntax:
 * precicef_read_sdata(
 *   INTEGER dataID,
 *   INTEGER valueIndex,
 *   DOUBLE PRECISION dataValue )
 *
 * IN:  dataID, valueIndex
 * OUT: dataValue
 *
 * @copydoc precice::SolverInterface::readScalarData
 *
 */
void precicef_read_sdata_(
    const int *dataID,
    const int *valueIndex,
    double *   dataValue);

/**
 * Fortran syntax:
 * precicef_map_write_data_from( INTEGER meshID )
 *
 * IN:  meshID
 * OUT: -
 *
 * @copydoc precice::SolverInterface::mapWriteDataFrom()
 *
 */
void precicef_map_write_data_from_(const int *meshID);

/**
 * Fortran syntax:
 * precicef_map_read_data_to( INTEGER meshID )
 *
 * IN:  meshID
 * OUT: -
 *
 * @copydoc precice::SolverInterface::mapReadDataTo()
 *
 */
void precicef_map_read_data_to_(const int *meshID);

/**
 * @brief Name of action for writing iteration checkpoint.
 *
 * Fortran syntax:
 * precicef_action_write_iter_checkpoint( CHARACTER nameAction(*) )
 */
void precicef_action_write_iter_checkp_(
    char *nameAction,
    int   lengthNameAction);

/**
 * @brief Name of action for writing initial data.
 *
 * FortranSyntax:
 * precicef_action_write_initial_data( CHARACTER nameAction(*) )
 */
void precicef_action_write_initial_data_(
    char *nameAction,
    int   lengthNameAction);

/**
 * @brief Name of action for reading iteration checkpoint.
 *
 * Fortran syntax:
 * precicef_action_read_iter_checkpoint( CHARACTER nameAction(*) )
 */
void precicef_action_read_iter_checkp_(
    char *nameAction,
    int   lengthNameAction);

void precicef_get_version_information_(
    char *versionInfo,
    int   lengthVersionInfo);

/** @name Experimental Data Access
 * These API functions are \b experimental and may change in future versions.
 */
///@{

/**
 * Fortran syntax:
 * precicef_set_mesh_access_region_(
 *   INTEGER          meshID,
 *   DOUBLE PRECISION bounding_box(dim*2))
 *
 * IN:  meshID, bounding_box
 * OUT: -
 *
 * @copydoc precice::SolverInterface::setMeshAccessRegion()
 */
void precicef_set_mesh_access_region_(
    const int     meshID,
    const double *boundingBox);

/**
 * Fortran syntax:
 * precicef_get_mesh_vertices_and_IDs_(
 *   INTEGER          meshID,
 *   INTEGER          size,
 *   INTEGER          ids(size),
 *   DOUBLE PRECISION coordinates(dim*size))
 *
 * IN:  meshID, size
 * OUT: ids, coordinates
 *
 * @copydoc precice::SolverInterface::getMeshVerticesAndIDs()
 */
void precicef_get_mesh_vertices_and_IDs_(
    const int meshID,
    const int size,
    int *     ids,
    double *  coordinates);

///@}

#ifdef __cplusplus
}
#endif
