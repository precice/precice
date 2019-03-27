#pragma once

#include "cplscheme/SharedPointer.hpp"
#include "com/SharedPointer.hpp"
#include "logging/Logger.hpp"
#include <set>
#include <list>
#include <Eigen/Core>

namespace precice {
  namespace impl {
    class SolverInterfaceImpl;
  }
}

namespace precice {
namespace impl {

/// Takes requests from clients and handles requests on server side.
class RequestManager
{
public:

  RequestManager (
    SolverInterfaceImpl&  solverInterfaceImpl,
    com::PtrCommunication clientServerCommunication,
    cplscheme::PtrCouplingScheme couplingScheme);

  /// Redirects all requests from client to corresponding handle methods.
  void handleRequests();

  /// Pings server.
  void requestPing();

  /// Requests initialization from server.
  void requestInitialize();

  /// Requests initialization of data from server.
  void requestInitialzeData();

  /// Requests advance from server.
  void requestAdvance ( double dt );

  /// Requests finalize from server.
  void requestFinalize();

  /// Requests fulfilled action from server.
  void requestFulfilledAction ( const std::string& action );

  /// Requests set position of solver mesh from server.
  int requestSetMeshVertex (
    int               meshID,
    Eigen::VectorXd&  position );

  /// Requests get size of vertices of preCICE mesh.
  int requestGetMeshVertexSize(int meshID);

  /// Requests reset of a preCICE mesh.
  void requestResetMesh(int meshID);

  /// Requests set vertex positions from server.
  void requestSetMeshVertices (
    int     meshID,
    int     size,
    double* positions,
    int*    ids );

  /// Requests get vertex positions from server.
  void requestGetMeshVertices (
    int     meshID,
    int     size,
    int*    ids,
    double* positions );

  /// Requests get vertex ids from server.
  void requestGetMeshVertexIDsFromPositions (
    int     meshID,
    int     size,
    double* positions,
    int*    ids );

  /// Requests set mesh edge from server.
  int requestSetMeshEdge (
    int meshID,
    int firstVertexID,
    int secondVertexID );

  /// Requests set mesh triangle from server.
  void requestSetMeshTriangle (
    int meshID,
    int firstEdgeID,
    int secondEdgeID,
    int thirdEdgeID );

  /// Requests set mesh triangle with edges from server.
  void requestSetMeshTriangleWithEdges (
    int meshID,
    int firstVertexID,
    int secondVertexID,
    int thirdVertexID );

  /// Requests set mesh quad from server.
  void requestSetMeshQuad (
    int meshID,
    int firstEdgeID,
    int secondEdgeID,
    int thirdEdgeID,
    int fourthEdgeID );

  /// Requests set mesh quad with edges from server.
  void requestSetMeshQuadWithEdges (
    int meshID,
    int firstVertexID,
    int secondVertexID,
    int thirdVertexID,
    int fourthVertexID );

  /// Requests write block scalar data from server.
  void requestWriteBlockScalarData (
    int           dataID,
    int           size,
    int*          valueIndices,
    const double* values );

  /// Requests write scalar data from server.
  void requestWriteScalarData (
    int           dataID,
    int           valueIndex,
    const double& value );

  /// Requests write block vector data from server.
  void requestWriteBlockVectorData (
    int           dataID,
    int           size,
    int*          valueIndices,
    const double* values );

  /// Requests write vector data from server.
  void requestWriteVectorData (
    int           dataID,
    int           valueIndex,
    const double* value );

  /// Requests read block scalar data from server.
  void requestReadBlockScalarData (
    int     dataID,
    int     size,
    int*    valueIndices,
    double* values );

  /// Requests read scalar data from server.
  void requestReadScalarData (
    int     dataID,
    int     valueIndex,
    double& value );

  /// Requests read block vector data from server.
  void requestReadBlockVectorData (
    int     dataID,
    int     size,
    int*    valueIndices,
    double* values );

  /// Requests read vector data from server.
  void requestReadVectorData (
    int     dataID,
    int     valueIndex,
    double* value );

  /// Requests write mapping data from server.
  void requestMapWriteDataFrom ( int fromMeshID );

  /// Requests read mapping data from server.
  void requestMapReadDataTo( int toMeshID );

private:

  /// IDs for requests from clients.
  enum Request {
    REQUEST_INITIALIZE,
    REQUEST_INITIALIZE_DATA,
    REQUEST_ADVANCE,
    REQUEST_FINALIZE,
    REQUEST_FULFILLED_ACTION,
    REQUEST_SET_MESH_VERTEX,
    REQUEST_GET_MESH_VERTEX_SIZE,
    REQUEST_RESET_MESH,
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
    REQUEST_PING // Used in tests only
  };

  logging::Logger _log{"impl::RequestManager"};

  SolverInterfaceImpl& _interface;

  com::PtrCommunication _com;

  cplscheme::PtrCouplingScheme _couplingScheme;

  /// Handles request initialize from client.
  void handleRequestInitialze ( const std::list<int>& clientRanks );

  /// Handles request initialize data from client.
  void handleRequestInitialzeData ( const std::list<int>& clientRanks );

  /// Handles request advance from client.
  void handleRequestAdvance ( const std::list<int>& clientRanks );

  /// Handles request finalize from client.
  void handleRequestFinalize();

  /// Handles request fulfilled action from client.
  void handleRequestFulfilledAction ( int rankSender );

  /// Handles request set mesh vertex from client.
  void handleRequestSetMeshVertex ( int rankSender );

  /// Handles request get mesh vertex size from client.
  void handleRequestGetMeshVertexSize(int rankSender);

  /// Handles request reset mesh from client.
  void handleRequestResetMesh(int rankSender);

  /// Handles request set vertex positions from client.
  void handleRequestSetMeshVertices ( int rankSender );

  /// Handles request get vertex positions from client.
  void handleRequestGetMeshVertices ( int rankSender );

  /// Handles request get vertex IDs from client.
  void handleRequestGetMeshVertexIDsFromPositions ( int rankSender );

  /// Handles request set mesh edge from client.
  void handleRequestSetMeshEdge ( int rankSender );

  /// Handles request set mesh triangle from client.
  void handleRequestSetMeshTriangle ( int rankSender );

  /// Handles request set mesh triangle with edges from client.
  void handleRequestSetMeshTriangleWithEdges ( int rankSender );

  /// Handles request set mesh quad from client.
  void handleRequestSetMeshQuad ( int rankSender );

  /// Handles request set mesh quad with edges from client.
  void handleRequestSetMeshQuadWithEdges ( int rankSender );

  /// Handles request write block scalar data from client.
  void handleRequestWriteBlockScalarData ( int rankSender );

  /// Handles request write scalar data from client.
  void handleRequestWriteScalarData ( int rankSender );

  /// Handles request write block vector data from client.
  void handleRequestWriteBlockVectorData ( int rankSender );

  /// Handles request write vector data from client.
  void handleRequestWriteVectorData ( int rankSender );

  /// Handles request read block scalar data from client.
  void handleRequestReadBlockScalarData ( int rankSender );

  /// Handles request read scalar data from client.
  void handleRequestReadScalarData ( int rankSender );

  /// Handles request read block vector data from client.
  void handleRequestReadBlockVectorData ( int rankSender );

  /// Handles request read vector data from client.
  void handleRequestReadVectorData ( int rankSender );

  /// Handles request map written data from client.
  void handleRequestMapWriteDataFrom ( const std::list<int>& clientRanks );

  /// Handles request map read data from client.
  void handleRequestMapReadDataTo ( const std::list<int>& clientRanks );

};

}} // namespace precice, impl
