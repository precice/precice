// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_IMPL_REQUESTMANAGER_HPP_
#define PRECICE_IMPL_REQUESTMANAGER_HPP_

#include "precice/ClosestMesh.hpp"
#include "precice/VoxelPosition.hpp"
#include "com/Communication.hpp"
#include "cplscheme/SharedPointer.hpp"
#include "tarch/logging/Log.h"
#include "utils/Dimensions.hpp"
#include <set>
#include <list>

namespace precice {
  namespace impl {
    class SolverInterfaceImpl;
  }
}

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace impl {

/**
 * @brief Takes requests from clients and handles requests on server side.
 */
class RequestManager
{
public:

  RequestManager (
    bool                  geometryMode,
    SolverInterfaceImpl&  solverInterfaceImpl,
    com::Communication::SharedPointer clientServerCommunication,
    cplscheme::PtrCouplingScheme couplingScheme);

  /**
   * @brief Redirects all requests from client to corresponding handle methods.
   */
  void handleRequests();

  /**
   * @brief Pings server.
   */
  void requestPing();

  /**
   * @brief Requests initialization from server.
   */
  void requestInitialize();

  /**
   * @brief Requests initialization of data from server.
   */
  void requestInitialzeData();

  /**
   * @brief Requests advance from server.
   */
  void requestAdvance ( double dt );

  /**
   * @brief Requests finalize from server.
   */
  void requestFinalize();

  /**
   * @brief Requests fulfilled action from server.
   */
  void requestFulfilledAction ( const std::string& action );

  /**
   * @brief Requests inquire position from server.
   */
  int requestInquirePosition (
    utils::DynVector&    point,
    const std::set<int>& meshIDs );

  /**
   * @brief Requests inquire closest mesh from server.
   */
  void requestInquireClosestMesh (
    utils::DynVector&    point,
    const std::set<int>& meshIDs,
    ClosestMesh&         closest );

  /**
   * @brief Requests inquire voxel position from server.
   */
  void requestInquireVoxelPosition (
    utils::DynVector&    voxelCenter,
    utils::DynVector&    voxelHalflengths,
    bool                 includeBoundaries,
    const std::set<int>& meshIDs,
    VoxelPosition&       voxelPosition );

  /**
   * @brief Requests set position of solver mesh from server.
   */
  int requestSetMeshVertex (
    int               meshID,
    utils::DynVector& position );

  /**
   * @brief Requests get size of vertices of preCICE mesh.
   */
  int requestGetMeshVertexSize(int meshID);

  /**
   * @brief Requests set vertex positions from server.
   */
  void requestSetMeshVertices (
    int     meshID,
    int     size,
    double* positions,
    int*    ids );

  /**
   * @brief Requests get vertex positions from server.
   */
  void requestGetMeshVertices (
    int     meshID,
    int     size,
    int*    ids,
    double* positions );

  /**
   * @brief Requests get vertex ids from server.
   */
  void requestGetMeshVertexIDsFromPositions (
    int     meshID,
    int     size,
    double* positions,
    int*    ids );

  /**
   * @brief Requests set mesh edge from server.
   */
  int requestSetMeshEdge (
    int meshID,
    int firstVertexID,
    int secondVertexID );

  /**
   * @brief Requests set mesh triangle from server.
   */
  void requestSetMeshTriangle (
    int meshID,
    int firstEdgeID,
    int secondEdgeID,
    int thirdEdgeID );

  /**
   * @brief Requests set mesh triangle with edges from server.
   */
  void requestSetMeshTriangleWithEdges (
    int meshID,
    int firstVertexID,
    int secondVertexID,
    int thirdVertexID );

  /**
   * @brief Requests set mesh quad from server.
   */
  void requestSetMeshQuad (
    int meshID,
    int firstEdgeID,
    int secondEdgeID,
    int thirdEdgeID,
    int fourthEdgeID );

  /**
   * @brief Requests set mesh quad with edges from server.
   */
  void requestSetMeshQuadWithEdges (
    int meshID,
    int firstVertexID,
    int secondVertexID,
    int thirdVertexID,
    int fourthVertexID );

  /**
   * @brief Requests write block scalar data from server.
   */
  void requestWriteBlockScalarData (
    int     dataID,
    int     size,
    int*    valueIndices,
    double* values );

  /**
   * @brief Requests write scalar data from server.
   */
  void requestWriteScalarData (
    int    dataID,
    int    valueIndex,
    double value );

  /**
   * @brief Requests write block vector data from server.
   */
  void requestWriteBlockVectorData (
    int     dataID,
    int     size,
    int*    valueIndices,
    double* values );

  /**
   * @brief Requests write vector data from server.
   */
  void requestWriteVectorData (
    int     dataID,
    int     valueIndex,
    double* value );

  /**
   * @brief Requests read block scalar data from server.
   */
  void requestReadBlockScalarData (
    int     dataID,
    int     size,
    int*    valueIndices,
    double* values );

  /**
   * @brief Requests read scalar data from server.
   */
  void requestReadScalarData (
    int     dataID,
    int     valueIndex,
    double& value );

  /**
   * @brief Requests read block vector data from server.
   */
  void requestReadBlockVectorData (
    int     dataID,
    int     size,
    int*    valueIndices,
    double* values );

  /**
   * @brief Requests read vector data from server.
   */
  void requestReadVectorData (
    int     dataID,
    int     valueIndex,
    double* value );

  /**
   * @brief Requests write mapping data from server.
   */
  void requestMapWriteDataFrom ( int fromMeshID );

  /**
   * @brief Requests read mapping data from server.
   */
  void requestMapReadDataTo( int toMeshID );

  /**
   * @brief Requests export mesh from server.
   */
  void requestExportMesh (
    const std::string& filenameSuffix,
    int                exportType );

private:

  /**
   * @brief IDs for requests from clients.
   */
  enum Request {
    REQUEST_INITIALIZE,
    REQUEST_INITIALIZE_DATA,
    REQUEST_ADVANCE,
    REQUEST_FINALIZE,
    REQUEST_FULFILLED_ACTION,
    REQUEST_INQUIRE_POSITION,
    REQUEST_INQUIRE_CLOSEST_MESH,
    REQUEST_INQUIRE_VOXEL_POSITION,
    REQUEST_SET_MESH_VERTEX,
    REQUEST_GET_MESH_VERTEX_SIZE,
    REQUEST_SET_MESH_VERTICES,
    REQUEST_GET_MESH_VERTICES,
    REQUEST_GET_MESH_VERTEX_IDS_FROM_POSITIONS,
    REQUEST_SET_MESH_EDGE,
    REQUEST_SET_MESH_TRIANGLE,
    REQUEST_SET_MESH_TRIANGLE_WITH_EDGES,
    REQUEST_SET_MESH_QUAD,
    REQUEST_SET_MESH_QUAD_WITH_EDGES,
    REQUEST_WRITE_SCALAR_DATA,
    REQUEST_WRITE_BLOCK_SCALAR_DATA,
    REQUEST_WRITE_VECTOR_DATA,
    REQUEST_WRITE_BLOCK_VECTOR_DATA,
    REQUEST_READ_SCALAR_DATA,
    REQUEST_READ_BLOCK_SCALAR_DATA,
    REQUEST_READ_VETOR_DATA,
    REQUEST_READ_BLOCK_VECTOR_DATA,
    REQUEST_MAP_WRITE_DATA_FROM,
    REQUEST_MAP_READ_DATA_TO,
    REQUEST_EXPORT_MESH,
    REQUEST_PING // Used in tests only
  };

  static tarch::logging::Log _log;

  bool _isGeometryMode;

  SolverInterfaceImpl& _interface;

  com::Communication::SharedPointer _com;

  cplscheme::PtrCouplingScheme _couplingScheme;

  // @brief Locks the next receive operation of the server to a specific client.
  int _lockServerToClient;

  /**
   * @brief Handles request initialize from client.
   */
  void handleRequestInitialze ( const std::list<int>& clientRanks );

  /**
   * @brief Handles request initialize data from client.
   */
  void handleRequestInitialzeData ( const std::list<int>& clientRanks );

  /**
   * @brief Handles request advance from client.
   */
  void handleRequestAdvance ( const std::list<int>& clientRanks );

  /**
   * @brief Handles request finalize from client.
   */
  void handleRequestFinalize();

  /**
   * @brief Handles request fulfilled action from client.
   */
  void handleRequestFulfilledAction ( int rankSender );

  /**
   * @brief Handles request inquire position from client.
   */
  void handleRequestInquirePosition ( int rankSender );

  /**
   * @brief Handles request inquire closest mesh from client.
   */
  void handleRequestInquireClosestMesh ( int rankSender );

  /**
   * @brief Handles request inquire voxel position from client.
   */
  void handleRequestInquireVoxelPosition ( int rankSender );

  /**
   * @brief Handles request set mesh vertex from client.
   */
  void handleRequestSetMeshVertex ( int rankSender );

  /**
   * @brief Handles request get mesh vertex size from client.
   */
  void handleRequestGetMeshVertexSize(int rankSender);

  /**
   * @brief Handles request set vertex positions from client.
   */
  void handleRequestSetMeshVertices ( int rankSender );

  /**
   * @brief Handles request get vertex positions from client.
   */
  void handleRequestGetMeshVertices ( int rankSender );

  /**
   * @brief Handles request get vertex IDs from client.
   */
  void handleRequestGetMeshVertexIDsFromPositions ( int rankSender );

  /**
   * @brief Handles request set mesh edge from client.
   */
  void handleRequestSetMeshEdge ( int rankSender );

  /**
   * @brief Handles request set mesh triangle from client.
   */
  void handleRequestSetMeshTriangle ( int rankSender );

  /**
   * @brief Handles request set mesh triangle with edges from client.
   */
  void handleRequestSetMeshTriangleWithEdges ( int rankSender );

  /**
   * @brief Handles request set mesh quad from client.
   */
  void handleRequestSetMeshQuad ( int rankSender );

  /**
   * @brief Handles request set mesh quad with edges from client.
   */
  void handleRequestSetMeshQuadWithEdges ( int rankSender );

  /**
   * @brief Handles request write block scalar data from client.
   */
  void handleRequestWriteBlockScalarData ( int rankSender );

  /**
   * @brief Handles request write scalar data from client.
   */
  void handleRequestWriteScalarData ( int rankSender );

  /**
   * @brief Handles request write block vector data from client.
   */
  void handleRequestWriteBlockVectorData ( int rankSender );

  /**
   * @brief Handles request write vector data from client.
   */
  void handleRequestWriteVectorData ( int rankSender );

  /**
   * @brief Handles request read block scalar data from client.
   */
  void handleRequestReadBlockScalarData ( int rankSender );

  /**
   * @brief Handles request read scalar data from client.
   */
  void handleRequestReadScalarData ( int rankSender );

  /**
   * @brief Handles request read block vector data from client.
   */
  void handleRequestReadBlockVectorData ( int rankSender );

  /**
   * @brief Handles request read vector data from client.
   */
  void handleRequestReadVectorData ( int rankSender );

  /**
   * @brief Handles request map written data from client.
   */
  void handleRequestMapWriteDataFrom ( const std::list<int>& clientRanks );

  /**
   * @brief Handles request map read data from client.
   */
  void handleRequestMapReadDataTo ( const std::list<int>& clientRanks );

  /**
   * @brief Handles request export mesh from client.
   */
  void handleRequestExportMesh ( int rankSender );

};

}} // namespace precice, impl

#endif /* PRECICE_IMPL_REQUESTMANAGER_HPP_ */
