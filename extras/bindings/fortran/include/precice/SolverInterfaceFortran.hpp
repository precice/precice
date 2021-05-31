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
 * @brief See precice::SolverInterface::SolverInterface() and ...::configure().
 *
 * Fortran syntax:
 * precicef_create(
 *   CHARACTER participantName(*),
 *   CHARACTER configFileName(*),
 *   INTEGER   solverProcessIndex,
 *   INTEGER   solverProcessSize )
 *
 * IN:  participantName, configFileName, solverProcessIndex, solverProcessSize
 * OUT: -
 */
void precicef_create_(
    const char *participantName,
    const char *configFileName,
    const int * solverProcessIndex,
    const int * solverProcessSize,
    int         lengthAccessorName,
    int         lengthConfigFileName);

/**
 * @brief See precice::SolverInterface::initialize().
 *
 * Fortran syntax:
 * precicef_initialize( DOUBLE PRECISION timstepLengthLimit )
 *
 * IN:  -
 * OUT: timestepLengthLimit
 */
void precicef_initialize_(double *timestepLengthLimit);

/**
 * @brief See precice::SolverInterface::initializeData().
 *
 * Fortran syntax:
 * precicef_intialize_data()
 *
 * IN: -
 * OUT: -
 */
void precicef_initialize_data_();

/**
 * @brief See precice::SolverInterface::advance().
 *
 * Fortran syntax:
 * precicef_advance( DOUBLE PRECISION timstepLengthLimit )
 *
 * IN:  timestepLengthLimit
 * OUT: timestepLengthLimit
 */
void precicef_advance_(double *timestepLengthLimit);

/**
 * @brief See precice::SolverInterface::finalize().
 *
 * Fortran syntax:
 * precicef_finalize();
 */
void precicef_finalize_();

/**
 * @brief See precice::SolverInterface::getDimensions().
 *
 * Fortran syntax:
 * precicef_get_dims( INTEGER dimensions )
 *
 * IN:  -
 * OUT: dimensions
 */
void precicef_get_dims_(int *dimensions);

/**
 * @brief See precice::SolverInterface::isOngoing().
 *
 * @deprecated Forwards to precicef_is_coupling_ongoing_
 *
 * Fortran syntax:
 * precicef_ongoing( INTEGER isOngoing )
 *
 * IN:  -
 * OUT: isOngoing(1:true, 0:false)
 */
[[deprecated("Use precicef_is_coupling_ongoing_() instead.")]] void precicef_ongoing_(int *isOngoing);

/**
 * @brief See precice::SolverInterface::isCouplingOngoing().
 *
 * Fortran syntax:
 * precicef_is_coupling_ongoing( INTEGER isOngoing )
 *
 * IN:  -
 * OUT: isOngoing(1:true, 0:false)
 */
void precicef_is_coupling_ongoing_(int *isOngoing);

/**
 * @brief See precice::SolverInterface::isWriteDataRequired().
 *
 * @deprecated Forwards to precicef_is_write_data_required_
 *
 * Fortran syntax:
 * precicef_write_data_required(
 *  DOUBLE PRECISION computedTimestepLength,
 *  INTEGER          isRequired )
 *
 * IN:  computedTimestepLength
 * OUT: isRequired(1:true, 0:false)
 */
[[deprecated("Use precicef_is_write_data_required_(...) with the same arguments instead.")]] void precicef_write_data_required_(
    const double *computedTimestepLength,
    int *         isRequired);

/**
 * @brief See precice::SolverInterface::isWriteDataRequired().
 *
 * Fortran syntax:
 * precicef_is_write_data_required(
 *  DOUBLE PRECISION computedTimestepLength,
 *  INTEGER          isRequired )
 *
 * IN:  computedTimestepLength
 * OUT: isRequired(1:true, 0:false)
 */
void precicef_is_write_data_required_(
    const double *computedTimestepLength,
    int *         isRequired);

/**
 * @brief See precice::SolverInterface::isReadDataAvailable().
 *
 * @deprecated Forwards to precicef_is_read_data_available_
 *
 * Fortran syntax:
 * precicef_read_data_available( INTEGER isAvailable );
 *
 * IN:  -
 * OUT: isAvailable(1:true, 0:false)
 */
[[deprecated("Use precicef_is_read_data_available_() instead.")]] void precicef_read_data_available_(int *isAvailable);

/**
 * @brief See precice::SolverInterface::isReadDataAvailable().
 *
 * Fortran syntax:
 * precicef_is_read_data_available( INTEGER isAvailable );
 *
 * IN:  -
 * OUT: isAvailable(1:true, 0:false)
 */
void precicef_is_read_data_available_(int *isAvailable);

/**
 * @brief See precice::SolverInterface::isTimeWindowComplete().
 *
 * Fortran syntax:
 * precicef_is_time_window_complete( INTEGER isComplete );
 *
 * IN:  -
 * OUT: isComplete(1:true, 0:false)
 */
void precicef_is_time_window_complete_(int *isComplete);

/**
 * @brief See precice::SolverInterface::hasToEvaluateSurrogateModel().
 *
 * Fortran syntax:
 * precicef_has_to_evaluate_surrogate_model( INTEGER hasToEvaluate );
 *
 * IN:  -
 * OUT: hasToEvaluate(1:true, 0:false)
 */
void precicef_has_to_evaluate_surrogate_model_(int *hasToEvaluate);

/**
 * @brief See precice::SolverInterface::hasToEvaluateFineModel().
 *
 * Fortran syntax:
 * precicef_has_to_evaluate_fine_model( INTEGER hasToEvaluate );
 *
 * IN:  -
 * OUT: hasToEvaluate(1:true, 0:false)
 */
void precicef_has_to_evaluate_fine_model_(int *hasToEvaluate);

/**
 * @brief See precice::SolverInterface::isActionRequired().
 *
 * @deprecated Forwards to precicef_is_action_required_
 *
 * Fortran syntax:
 * precicef_action_required(
 *   CHARACTER action(*),
 *   INTEGER   isRequired )
 *
 * IN:  action
 * OUT: isRequired(1:true, 0:false)
 */
[[deprecated("Use precicef_is_action_required_(...) with the same arguments instead.")]] void precicef_action_required_(
    const char *action,
    int *       isRequired,
    int         lengthAction);

/**
 * @brief See precice::SolverInterface::isActionRequired().
 *
 * Fortran syntax:
 * precicef_is_action_required(
 *   CHARACTER action(*),
 *   INTEGER   isRequired )
 *
 * IN:  action
 * OUT: isRequired(1:true, 0:false)
 */
void precicef_is_action_required_(
    const char *action,
    int *       isRequired,
    int         lengthAction);

/**
 * @brief See precice::SolverInterface::markActionFulfilled().
 *
 * Fortran syntax:
 * precicef_mark_action_fulfilled( CHARACTER action(*) )
 *
 * IN:  action
 * OUT: -
 */
void precicef_mark_action_fulfilled_(
    const char *action,
    int         lengthAction);

/**
 * @brief See precice::SolverInterface::hasMesh().
 *
 * Fortran syntax:
 * precicef_has_mesh(
 *   CHARACTER meshName(*),
 *   INTEGER   hasMesh )
 *
 * IN:  meshName
 * OUT: hasMesh(1:true, 0:false)
 */
void precicef_has_mesh_(
    const char *meshName,
    int *       hasMesh,
    int         lengthMeshName);

/**
 * @brief See precice::SolverInterface::getMeshID().
 *
 * Fortran syntax:
 * precicef_get_mesh_id(
 *   CHARACTER meshName(*),
 *   INTEGER   meshID )
 *
 * IN:  meshName
 * OUT: meshID
 */
void precicef_get_mesh_id_(
    const char *meshName,
    int *       meshID,
    int         lengthMeshName);

/**
 * @brief See precice::SolverInterface::hasData().
 *
 * Fortran syntax:
 * precicef_has_data(
 *   CHARACTER dataName(*),
 *   INTEGER   meshID,
 *   INTEGER   hasData)
 *
 * IN:  dataName
 * IN:  meshID
 * OUT: hasData(1:true, 0:false)
 */
void precicef_has_data_(
    const char *dataName,
    const int * meshID,
    int *       hasData,
    int         lengthDataName);

/**
 * @brief See precice::SolverInterface::getDataID().
 *
 * The given name (dataName) has to be one of the names specified in the
 * configuration file. The data id obtained can be used to read and write
 * data to and from the coupling mesh.
 *
 * Fortran syntax:
 * precicef_get_data_id(
 *   CHARACTER dataName(*),
 *   INTEGER   meshID,
 *   INTEGER   dataID,
)
 *
 * IN:  dataName
 * IN:  meshID
 * OUT: dataID
 */
void precicef_get_data_id_(
    const char *dataName,
    const int * meshID,
    int *       dataID,
    int         lengthDataName);

/**
 * @brief See precice::SolverInterface::setMeshVertex().
 *
 * Fortran syntax:
 * precicef_set_vertex(
 *   INTEGER          meshID,
 *   DOUBLE PRECISION position(dim),
 *   INTEGER          vertexID )
 *
 * IN:  meshID, position
 * OUT: vertexID
 */
void precicef_set_vertex_(
    const int *   meshID,
    const double *position,
    int *         vertexID);

/**
 * @brief See precice::SolverInterface::getMeshVertexSize().
 *
 * Fortran syntax:
 * precicef_get_mesh_vertex_size(
 *   INTEGER meshID,
 *   INTEGER meshSize )
 *
 * IN:  meshID
 * OUT: meshSize
 */
void precicef_get_mesh_vertex_size_(
    const int *meshID,
    int *      meshSize);

/**
 * @brief See precice::SolverInterface::setMeshVertices().
 *
 * Fortran syntax:
 * precicef_set_vertices(
 *   INTEGER          meshID,
 *   INTEGER          size,
 *   DOUBLE PRECISION positions(dim*size),
 *   INTEGER          positionIDs(size) )
 *
 * IN:  meshID, size, positions
 * OUT: positionIDs
 */
void precicef_set_vertices_(
    const int *meshID,
    const int *size,
    double *   positions,
    int *      positionIDs);

/**
 * @brief See precice::SolverInterface::getMeshVertices().
 *
 * Fortran syntax:
 * precicef_get_vertices(
 *   INTEGER          meshID,
 *   INTEGER          size,
 *   INTEGER          ids(size)
 *   DOUBLE PRECISION positions(dim*size))
 *
 * IN:  meshID, size, ids
 * OUT: positions
 */
void precicef_get_vertices_(
    const int *meshID,
    const int *size,
    int *      ids,
    double *   positions);

/**
 * @brief See precice::SolverInterface::getMeshVertexIDsFromPositions().
 *
 * Fortran syntax:
 * precicef_get_vertices(
 *   INTEGER          meshID,
 *   INTEGER          size,
 *   DOUBLE PRECISION positions(dim*size),
 *   INTEGER          ids(size))
 *
 * IN:  meshID, size, positions
 * OUT: ids
 */
void precicef_get_vertex_ids_from_positions_(
    const int *meshID,
    const int *size,
    double *   positions,
    int *      ids);

/**
 * @brief See precice::SolverInterface::setMeshEdge().
 *
 * Fortran syntax:
 * precicef_set_edge(
 *   INTEGER meshID,
 *   INTEGER firstVertexID,
 *   INTEGER secondVertexID,
 *   INTEGER edgeID )
 *
 * IN:  meshID, firstVertexID, secondVertexID
 * OUT: edgeID
 */
void precicef_set_edge_(
    const int *meshID,
    const int *firstVertexID,
    const int *secondVertexID,
    int *      edgeID);

/**
 * @brief See precice::SolverInterface::setMeshTriangle().
 *
 * Fortran syntax:
 * precicef_set_triangle(
 *   INTEGER meshID,
 *   INTEGER firstEdgeID,
 *   INTEGER secondEdgeID,
 *   INTEGER thirdEdgeID )
 *
 * IN:  meshID, firstEdgeID, secondEdgeID, thirdEdgeID
 * OUT: -
 */
void precicef_set_triangle_(
    const int *meshID,
    const int *firstEdgeID,
    const int *secondEdgeID,
    const int *thirdEdgeID);

/**
 * @brief See precice::SolverInterface::setMeshTriangleWithEdges().
 *
 * Fortran syntax:
 * precicef_set_triangle_we(
 *   INTEGER meshID,
 *   INTEGER firstVertexID,
 *   INTEGER secondVertexID,
 *   INTEGER thirdVertexID )
 *
 * IN:  meshID, firstVertexID, secondVertexID, thirdVertexID
 * OUT: -
 */
void precicef_set_triangle_we_(
    const int *meshID,
    const int *firstVertexID,
    const int *secondVertexID,
    const int *thirdVertexID);

/**
 * @brief See precice::SolverInterface::setMeshQuad().
 *
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
 */
void precicef_set_quad_(
    const int *meshID,
    const int *firstEdgeID,
    const int *secondEdgeID,
    const int *thirdEdgeID,
    const int *fourthEdgeID);

/**
 * @brief See precice::SolverInterface::setMeshQuadWithEdges().
 *
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
 */
void precicef_set_quad_we_(
    const int *meshID,
    const int *firstVertexID,
    const int *secondVertexID,
    const int *thirdVertexID,
    const int *fourthVertexID);

/**
 * @brief See precice::SolverInterface::writeBlockVectorData.
 *
 * Fortran syntax:
 * precicef_write_bvdata(
 *   INTEGER dataID,
 *   INTEGER size,
 *   INTEGER valueIndices,
 *   DOUBLE PRECISION values(dim*size) )
 *
 * IN:  dataID, size, valueIndices, values
 * OUT: -
 */
void precicef_write_bvdata_(
    const int *dataID,
    const int *size,
    int *      valueIndices,
    double *   values);

/**
 * @brief precice::SolverInterface::writeVectorData.
 *
 * Fortran syntax:
 * precicef_write_vdata(
 *   INTEGER dataID,
 *   INTEGER valueIndex,
 *   DOUBLE PRECISION dataValue(dim) )
 *
 * IN:  dataID, valueIndex, dataValue
 * OUT: -
 */
void precicef_write_vdata_(
    const int *   dataID,
    const int *   valueIndex,
    const double *dataValue);

/**
 * @brief See precice::SolverInterface::writeBlockScalarData.
 *
 * Fortran syntax:
 * precicef_write_bsdata(
 *   INTEGER dataID,
 *   INTEGER size,
 *   INTEGER valueIndices,
 *   DOUBLE PRECISION values(size) )
 *
 * IN:  dataID, size, valueIndices, values
 * OUT: -
 */
void precicef_write_bsdata_(
    const int *dataID,
    const int *size,
    int *      valueIndices,
    double *   values);

/**
 * @brief precice::SolverInterface::writeScalarData.
 *
 * Fortran syntax:
 * precicef_write_sdata(
 *   INTEGER dataID,
 *   INTEGER valueIndex,
 *   DOUBLE PRECISION dataValue )
 *
 * IN:  dataID, valueIndex, dataValue
 * OUT: -
 */
void precicef_write_sdata_(
    const int *   dataID,
    const int *   valueIndex,
    const double *dataValue);

/**
 * @brief See precice::SolverInterface::readBlockVectorData.
 *
 * Fortran syntax:
 * precicef_read_bvdata(
 *   INTEGER dataID,
 *   INTEGER size,
 *   INTEGER valueIndices,
 *   DOUBLE PRECISION values(dim*size) )
 *
 * IN:  dataID, size, valueIndices
 * OUT: values
 */
void precicef_read_bvdata_(
    const int *dataID,
    const int *size,
    int *      valueIndices,
    double *   values);

/**
 * @brief precice::SolverInterface::readVectorData.
 *
 * Fortran syntax:
 * precicef_read_vdata(
 *   INTEGER dataID,
 *   INTEGER valueIndex,
 *   DOUBLE PRECISION dataValue(dim) )
 *
 * IN:  dataID, valueIndex
 * OUT: dataValue
 */
void precicef_read_vdata_(
    const int *dataID,
    const int *valueIndex,
    double *   dataValue);

/**
 * @brief See precice::SolverInterface::readBlockScalarData.
 *
 * Fortran syntax:
 * precicef_read_bsdata(
 *   INTEGER dataID,
 *   INTEGER size,
 *   INTEGER valueIndices,
 *   DOUBLE PRECISION values(size) )
 *
 * IN:  dataID, size, valueIndices
 * OUT: values
 */
void precicef_read_bsdata_(
    const int *dataID,
    const int *size,
    int *      valueIndices,
    double *   values);

/**
 * @brief precice::SolverInterface::readScalarData.
 *
 * Fortran syntax:
 * precicef_read_sdata(
 *   INTEGER dataID,
 *   INTEGER valueIndex,
 *   DOUBLE PRECISION dataValue )
 *
 * IN:  dataID, valueIndex
 * OUT: dataValue
 */
void precicef_read_sdata_(
    const int *dataID,
    const int *valueIndex,
    double *   dataValue);

/**
 * @brief See precice::SolverInterface::mapWriteDataFrom().
 *
 * Fortran syntax:
 * precicef_map_write_data_from( INTEGER meshID )
 *
 * IN:  meshID
 * OUT: -
 */
void precicef_map_write_data_from_(const int *meshID);

/**
 * @brief See precice::SolverInterface::mapReadDataTo().
 *
 * Fortran syntax:
 * precicef_map_read_data_to( INTEGER meshID )
 *
 * IN:  meshID
 * OUT: -
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

#ifdef __cplusplus
}
#endif
