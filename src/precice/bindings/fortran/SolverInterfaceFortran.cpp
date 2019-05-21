#include "SolverInterfaceFortran.hpp"
#include "precice/impl/SolverInterfaceImpl.hpp" // TODO remove this?
#include <iostream>
#include <string>
#include "logging/Logger.hpp"
#include "precice/SolverInterface.hpp"
#include "utils/assertion.hpp"

using namespace std;

static precice::impl::SolverInterfaceImpl* impl = nullptr;

static precice::logging::Logger _log ("SolverInterfaceFortran");

static std::string errormsg = "preCICE has not been created properly. Be sure to call \"precicef_create\" before any other call to preCICE.";

namespace precice {
  namespace impl {
    /**
     * @brief Returns length of string without trailing whitespaces.
     */
    int strippedLength( const char* string, int length );
  }
}

void precicef_create_
(
  const char* participantName,
  const char* configFileName,
  const int*  solverProcessIndex,
  const int*  solverProcessSize,
  int   lengthAccessorName,
  int   lengthConfigFileName )
{
  //cout << "lengthAccessorName: " << lengthAccessorName << '\n';
  //cout << "lengthConfigFileName: " << lengthConfigFileName << '\n';
  //cout << "solverProcessIndex: " << *solverProcessIndex << '\n';
  //cout << "solverProcessSize: " << *solverProcessSize << '\n';
  int strippedLength = precice::impl::strippedLength(participantName,lengthAccessorName);
  string stringAccessorName(participantName, strippedLength);
  strippedLength = precice::impl::strippedLength(configFileName,lengthConfigFileName);
  string stringConfigFileName(configFileName, strippedLength);
  //cout << "Accessor: " << stringAccessorName << "!" << '\n';
  //cout << "Config  : " << stringConfigFileName << "!" << '\n';
  impl = new precice::impl::SolverInterfaceImpl (stringAccessorName,
         *solverProcessIndex, *solverProcessSize, false);
  impl->configure(stringConfigFileName);
}

void precicef_initialize_
(
  double* timestepLengthLimit)
{
  CHECK(impl != nullptr,errormsg);
  *timestepLengthLimit = impl->initialize();
}

void precicef_initialize_data_()
{
  CHECK(impl != nullptr,errormsg);
  impl->initializeData();
}

void precicef_advance_
(
  double* timestepLengthLimit)
{
  CHECK(impl != nullptr,errormsg);
  *timestepLengthLimit = impl->advance(*timestepLengthLimit);
}

void precicef_finalize_()
{
  CHECK(impl != nullptr,errormsg);
  impl->finalize();
  delete impl;
}

void precicef_get_dims_
(
  int* dimensions )
{
  CHECK(impl != nullptr,errormsg);
  *dimensions = impl->getDimensions();
}

void precicef_ongoing_
(
  int* isOngoing )
{
  CHECK(impl != nullptr,errormsg);
  if (impl->isCouplingOngoing()){
    *isOngoing = 1;
  }
  else {
    *isOngoing = 0;
  }
}

void precicef_write_data_required_
(
  const double* computedTimestepLength,
  int*          isRequired )
{
  CHECK(impl != nullptr,errormsg);
  if (impl->isWriteDataRequired(*computedTimestepLength)){
    *isRequired = 1;
  }
  else {
    *isRequired = 0;
  }
}

void precicef_read_data_available_
(
  int* isAvailable )
{
  CHECK(impl != nullptr,errormsg);
  if (impl->isReadDataAvailable()){
    *isAvailable = 1;
  }
  else {
    *isAvailable = 0;
  }
}

void precicef_is_timestep_complete_(
    int* isComplete )
{
  CHECK(impl != nullptr,errormsg);
  if (impl->isTimestepComplete()){
    *isComplete = 1;
  }
  else {
    *isComplete = 0;
  }
}

void precicef_has_to_evaluate_surrogate_model_(
    int* hasToEvaluate)
{
  CHECK(impl != nullptr,errormsg);
  if (impl->hasToEvaluateSurrogateModel()){
    *hasToEvaluate = 1;
  }
  else {
    *hasToEvaluate = 0;
  }
}

void precicef_has_to_evaluate_fine_model_(
    int* hasToEvaluate)
{
  CHECK(impl != nullptr,errormsg);
  if (impl->hasToEvaluateFineModel()){
    *hasToEvaluate = 1;
  }
  else {
    *hasToEvaluate = 0;
  }
}

void precicef_action_required_
(
  const char* action,
  int*        isRequired,
  int         lengthAction )
{
  CHECK(impl != nullptr,errormsg);
  //assertion(lengthAction > 1);
  //std::cout << "lengthAction: " << lengthAction << '\n';
  //std::cout << "Action:";
  //for (int i=0; i < lengthAction; i++){
  //  std::cout << " a[" << i << "]=\"" << action[i] << "\"";
  //}
  //std::cout << '\n';
  int strippedLength = precice::impl::strippedLength(action, lengthAction);
  //std::cout << "strippedLength: " << strippedLength << '\n';
  //assertion(strippedLength > 1);
  string stringAction(action, strippedLength);
  if (impl->isActionRequired(stringAction)){
    *isRequired = 1;
  }
  else {
    *isRequired = 0;
  }
}

void precicef_fulfilled_action_
(
  const char* action,
  int         lengthAction )
{
  CHECK(impl != nullptr,errormsg);
  int strippedLength = precice::impl::strippedLength(action, lengthAction);
  string stringAction(action, strippedLength);
  impl->fulfilledAction(stringAction);
}

void precicef_has_mesh_(
  const char* meshName,
  int*        hasMesh,
  int         lengthMeshName)
{
  CHECK(impl != nullptr,errormsg);
  int strippedLength = precice::impl::strippedLength(meshName, lengthMeshName);
  string stringMeshName(meshName, strippedLength);
  if (impl->hasMesh(meshName)){
    *hasMesh = 1;
  }
  else {
    *hasMesh = 0;
  }
}

void precicef_get_mesh_id_
(
  const char* meshName,
  int*        meshID,
  int         lengthMeshName )
{
  CHECK(impl != nullptr,errormsg);
  int strippedLength = precice::impl::strippedLength(meshName, lengthMeshName);
  string stringMeshName(meshName, strippedLength);
  *meshID = impl->getMeshID(stringMeshName);
}

void precicef_has_data_
(
  const char* dataName,
  const int*  meshID,
  int*        hasData,
  int         lengthDataName)
{
  CHECK(impl != nullptr,errormsg);
  int strippedLength = precice::impl::strippedLength(dataName, lengthDataName);
  string stringDataName(dataName, strippedLength);
  if (impl->hasData(stringDataName, *meshID)){
    *hasData = 1;
  }
  else {
    *hasData = 0;
  }
}

void precicef_get_data_id_
(
  const char* dataName,
  const int*  meshID,
  int*        dataID,
  int         lengthDataName
)
{
  CHECK(impl != nullptr,errormsg);
  int strippedLength = precice::impl::strippedLength(dataName, lengthDataName);
  string stringDataName(dataName, strippedLength);
  *dataID = impl->getDataID(stringDataName, *meshID);
}

void precicef_set_vertex_
(
  const int*    meshID,
  const double* position,
  int*          vertexID )
{
  CHECK(impl != nullptr,errormsg);
  *vertexID = impl->setMeshVertex(*meshID, position);
}

void precicef_get_mesh_vertex_size_(
    const int* meshID,
    int* meshSize)
{
  CHECK(impl != nullptr,errormsg);
  *meshSize = impl->getMeshVertexSize(*meshID);
}

void precicef_set_vertices_
(
  const int*    meshID,
  const int*    size,
  double*       positions,
  int*          positionIDs )
{
  CHECK(impl != nullptr,errormsg);
  impl->setMeshVertices(*meshID, *size, positions, positionIDs);
}

void precicef_get_vertices_(
  const int* meshID,
  const int* size,
  int*       ids,
  double*    positions )
{
  CHECK(impl != nullptr,errormsg);
  impl->getMeshVertices(*meshID, *size, ids, positions);
}

void precicef_get_vertex_ids_from_positions_(
  const int* meshID,
  const int* size,
  double*    positions,
  int*       ids )
{
  CHECK(impl != nullptr,errormsg);
  impl->getMeshVertexIDsFromPositions(*meshID, *size, positions, ids);
}

void precicef_set_edge_
(
  const int* meshID,
  const int* firstVertexID,
  const int* secondVertexID,
  int* edgeID )
{
  CHECK(impl != nullptr,errormsg);
  *edgeID = impl->setMeshEdge(*meshID, *firstVertexID, *secondVertexID);
}

void precicef_set_triangle_
(
  const int* meshID,
  const int* firstEdgeID,
  const int* secondEdgeID,
  const int* thirdEdgeID )
{
  CHECK(impl != nullptr,errormsg);
  impl->setMeshTriangle(*meshID, *firstEdgeID, *secondEdgeID, *thirdEdgeID);
}

void precicef_set_triangle_we_
(
  const int* meshID,
  const int* firstVertexID,
  const int* secondVertexID,
  const int* thirdVertexID )
{
  CHECK(impl != nullptr,errormsg);
  impl->setMeshTriangleWithEdges(*meshID, *firstVertexID, *secondVertexID, *thirdVertexID);
}

void precicef_set_quad_(
  const int* meshID,
  const int* firstEdgeID,
  const int* secondEdgeID,
  const int* thirdEdgeID,
  const int* fourthEdgeID)
{
  CHECK(impl != nullptr,errormsg);
  impl->setMeshQuad(*meshID, *firstEdgeID, *secondEdgeID, *thirdEdgeID, *fourthEdgeID);
}

void precicef_set_quad_we_(
  const int* meshID,
  const int* firstVertexID,
  const int* secondVertexID,
  const int* thirdVertexID,
  const int* fourthVertexID)
{
  CHECK(impl != nullptr,errormsg);
  impl->setMeshQuadWithEdges(*meshID, *firstVertexID, *secondVertexID, *thirdVertexID, *fourthVertexID);
}

void precicef_write_bvdata_
(
  const int* dataID,
  const int* size,
  int*       valueIndices,
  double*    values )
{
  CHECK(impl != nullptr,errormsg);
  impl->writeBlockVectorData(*dataID, *size, valueIndices, values);
}

void precicef_write_vdata_
(
  const int*    dataID,
  const int*    valueIndex,
  const double* dataValue )
{
  CHECK(impl != nullptr,errormsg);
  impl->writeVectorData(*dataID, *valueIndex, dataValue);
}

void precicef_write_bsdata_
(
  const int* dataID,
  const int* size,
  int*       valueIndices,
  double*    values )
{
  CHECK(impl != nullptr,errormsg);
  impl->writeBlockScalarData(*dataID, *size, valueIndices, values);
}

void precicef_write_sdata_
(
  const int*    dataID,
  const int*    valueIndex,
  const double* dataValue )
{
  CHECK(impl != nullptr,errormsg);
  impl->writeScalarData(*dataID, *valueIndex, *dataValue);
}

void precicef_read_bvdata_
(
  const int* dataID,
  const int* size,
  int*       valueIndices,
  double*    values )
{
  CHECK(impl != nullptr,errormsg);
  impl->readBlockVectorData(*dataID, *size, valueIndices, values);
}

void precicef_read_vdata_
(
  const int* dataID,
  const int* valueIndex,
  double*    dataValue )
{
  CHECK(impl != nullptr,errormsg);
  impl->readVectorData(*dataID, *valueIndex, dataValue);
}

void precicef_read_bsdata_
(
  const int* dataID,
  const int* size,
  int*       valueIndices,
  double*    values )
{
  CHECK(impl != nullptr,errormsg);
  impl->readBlockScalarData(*dataID, *size, valueIndices, values);
}

void precicef_read_sdata_
(
  const int* dataID,
  const int* valueIndex,
  double*    dataValue )
{
  CHECK(impl != nullptr,errormsg);
  impl->readScalarData(*dataID, *valueIndex, *dataValue);
}

void precicef_map_write_data_from_
(
  const int* meshID )
{
  CHECK(impl != nullptr,errormsg);
  impl->mapWriteDataFrom(*meshID);
}

void precicef_map_read_data_to_
(
  const int* meshID )
{
  CHECK(impl != nullptr,errormsg);
  impl->mapReadDataTo(*meshID);
}

int precice::impl::strippedLength
(
  const char* string,
  int         length )
{
  int i=length-1;
  while (((string[i] == ' ') || (string[i] == 0)) && (i >= 0)){
    i--;
  }
  return i+1;
}

void precicef_action_write_iter_checkp_
(
  char*  nameAction,
  int lengthNameAction )
{
  const std::string& name = precice::constants::actionWriteIterationCheckpoint();
  assertion(name.size() < (size_t) lengthNameAction, name.size(), lengthNameAction);
  for (size_t i=0; i < name.size(); i++){
    nameAction[i] = name[i];
  }
}

void precicef_action_write_initial_data_(
  char*  nameAction,
  int lengthNameAction )
{
  const std::string& name = precice::constants::actionWriteInitialData();
  assertion(name.size() < (size_t) lengthNameAction, name.size(), lengthNameAction);
  for (size_t i=0; i < name.size(); i++){
    nameAction[i] = name[i];
  }
}

void precicef_action_read_iter_checkp_
(
  char*  nameAction,
  int lengthNameAction )
{
  const std::string& name = precice::constants::actionReadIterationCheckpoint();
  assertion(name.size() < (size_t) lengthNameAction, name.size(), lengthNameAction);
  for (size_t i=0; i < name.size(); i++){
    nameAction[i] = name[i];
  }
}

