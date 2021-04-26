#include "precice/SolverInterfaceFortran.hpp"
#include <iostream>
#include <memory>
#include <stddef.h>
#include <string>
#include "logging/LogMacros.hpp"
#include "logging/Logger.hpp"
#include "precice/SolverInterface.hpp"
#include "precice/impl/versions.hpp"
#include "utils/assertion.hpp"

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
#endif

using namespace std;

static std::unique_ptr<precice::SolverInterface> impl = nullptr;

static precice::logging::Logger _log("SolverInterfaceFortran");

static std::string errormsg = "preCICE has not been created properly. Be sure to call \"precicef_create\" before any other call to preCICE.";

namespace precice {
namespace impl {
/**
     * @brief Returns length of string without trailing whitespaces.
     */
int strippedLength(const char *string, int length);
} // namespace impl
} // namespace precice

void precicef_create_(
    const char *participantName,
    const char *configFileName,
    const int * solverProcessIndex,
    const int * solverProcessSize,
    int         lengthAccessorName,
    int         lengthConfigFileName)
{
  //cout << "lengthAccessorName: " << lengthAccessorName << '\n';
  //cout << "lengthConfigFileName: " << lengthConfigFileName << '\n';
  //cout << "solverProcessIndex: " << *solverProcessIndex << '\n';
  //cout << "solverProcessSize: " << *solverProcessSize << '\n';
  int    strippedLength = precice::impl::strippedLength(participantName, lengthAccessorName);
  string stringAccessorName(participantName, strippedLength);
  strippedLength = precice::impl::strippedLength(configFileName, lengthConfigFileName);
  string stringConfigFileName(configFileName, strippedLength);
  //cout << "Accessor: " << stringAccessorName << "!" << '\n';
  //cout << "Config  : " << stringConfigFileName << "!" << '\n';
  impl.reset(new precice::SolverInterface(stringAccessorName,
                                          stringConfigFileName,
                                          *solverProcessIndex, *solverProcessSize));
}

void precicef_initialize_(
    double *timestepLengthLimit)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  *timestepLengthLimit = impl->initialize();
}

void precicef_initialize_data_()
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->initializeData();
}

void precicef_advance_(
    double *timestepLengthLimit)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  *timestepLengthLimit = impl->advance(*timestepLengthLimit);
}

void precicef_finalize_()
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->finalize();
  impl.reset();
}

void precicef_get_dims_(
    int *dimensions)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  *dimensions = impl->getDimensions();
}

void precicef_ongoing_(
    int *isOngoing)
{
  precicef_is_coupling_ongoing_(isOngoing);
}

void precicef_is_coupling_ongoing_(
    int *isOngoing)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  if (impl->isCouplingOngoing()) {
    *isOngoing = 1;
  } else {
    *isOngoing = 0;
  }
}

void precicef_write_data_required_(
    const double *computedTimestepLength,
    int *         isRequired)
{
  precicef_is_write_data_required_(computedTimestepLength, isRequired);
}

void precicef_is_write_data_required_(
    const double *computedTimestepLength,
    int *         isRequired)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  if (impl->isWriteDataRequired(*computedTimestepLength)) {
    *isRequired = 1;
  } else {
    *isRequired = 0;
  }
}

void precicef_read_data_available_(
    int *isAvailable)
{
  precicef_is_read_data_available_(isAvailable);
}

void precicef_is_read_data_available_(
    int *isAvailable)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  if (impl->isReadDataAvailable()) {
    *isAvailable = 1;
  } else {
    *isAvailable = 0;
  }
}

void precicef_is_time_window_complete_(
    int *isComplete)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  if (impl->isTimeWindowComplete()) {
    *isComplete = 1;
  } else {
    *isComplete = 0;
  }
}

void precicef_has_to_evaluate_surrogate_model_(
    int *hasToEvaluate)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  if (impl->hasToEvaluateSurrogateModel()) {
    *hasToEvaluate = 1;
  } else {
    *hasToEvaluate = 0;
  }
}

void precicef_has_to_evaluate_fine_model_(
    int *hasToEvaluate)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  if (impl->hasToEvaluateFineModel()) {
    *hasToEvaluate = 1;
  } else {
    *hasToEvaluate = 0;
  }
}

void precicef_action_required_(
    const char *action,
    int *       isRequired,
    int         lengthAction)
{
  precicef_is_action_required_(action, isRequired, lengthAction);
}

void precicef_is_action_required_(
    const char *action,
    int *       isRequired,
    int         lengthAction)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  //PRECICE_ASSERT(lengthAction > 1);
  //std::cout << "lengthAction: " << lengthAction << '\n';
  //std::cout << "Action:";
  //for (int i=0; i < lengthAction; i++){
  //  std::cout << " a[" << i << "]=\"" << action[i] << "\"";
  //}
  //std::cout << '\n';
  int strippedLength = precice::impl::strippedLength(action, lengthAction);
  //std::cout << "strippedLength: " << strippedLength << '\n';
  //PRECICE_ASSERT(strippedLength > 1);
  string stringAction(action, strippedLength);
  if (impl->isActionRequired(stringAction)) {
    *isRequired = 1;
  } else {
    *isRequired = 0;
  }
}

void precicef_mark_action_fulfilled_(
    const char *action,
    int         lengthAction)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  int    strippedLength = precice::impl::strippedLength(action, lengthAction);
  string stringAction(action, strippedLength);
  impl->markActionFulfilled(stringAction);
}

void precicef_has_mesh_(
    const char *meshName,
    int *       hasMesh,
    int         lengthMeshName)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  int    strippedLength = precice::impl::strippedLength(meshName, lengthMeshName);
  string stringMeshName(meshName, strippedLength);
  if (impl->hasMesh(meshName)) {
    *hasMesh = 1;
  } else {
    *hasMesh = 0;
  }
}

void precicef_get_mesh_id_(
    const char *meshName,
    int *       meshID,
    int         lengthMeshName)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  int    strippedLength = precice::impl::strippedLength(meshName, lengthMeshName);
  string stringMeshName(meshName, strippedLength);
  *meshID = impl->getMeshID(stringMeshName);
}

void precicef_has_data_(
    const char *dataName,
    const int * meshID,
    int *       hasData,
    int         lengthDataName)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  int    strippedLength = precice::impl::strippedLength(dataName, lengthDataName);
  string stringDataName(dataName, strippedLength);
  if (impl->hasData(stringDataName, *meshID)) {
    *hasData = 1;
  } else {
    *hasData = 0;
  }
}

void precicef_get_data_id_(
    const char *dataName,
    const int * meshID,
    int *       dataID,
    int         lengthDataName)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  int    strippedLength = precice::impl::strippedLength(dataName, lengthDataName);
  string stringDataName(dataName, strippedLength);
  *dataID = impl->getDataID(stringDataName, *meshID);
}

void precicef_set_vertex_(
    const int *   meshID,
    const double *position,
    int *         vertexID)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  *vertexID = impl->setMeshVertex(*meshID, position);
}

void precicef_get_mesh_vertex_size_(
    const int *meshID,
    int *      meshSize)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  *meshSize = impl->getMeshVertexSize(*meshID);
}

void precicef_set_vertices_(
    const int *meshID,
    const int *size,
    double *   positions,
    int *      positionIDs)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->setMeshVertices(*meshID, *size, positions, positionIDs);
}

void precicef_get_vertices_(
    const int *meshID,
    const int *size,
    int *      ids,
    double *   positions)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->getMeshVertices(*meshID, *size, ids, positions);
}

void precicef_get_vertex_ids_from_positions_(
    const int *meshID,
    const int *size,
    double *   positions,
    int *      ids)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->getMeshVertexIDsFromPositions(*meshID, *size, positions, ids);
}

void precicef_set_edge_(
    const int *meshID,
    const int *firstVertexID,
    const int *secondVertexID,
    int *      edgeID)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  *edgeID = impl->setMeshEdge(*meshID, *firstVertexID, *secondVertexID);
}

void precicef_set_triangle_(
    const int *meshID,
    const int *firstEdgeID,
    const int *secondEdgeID,
    const int *thirdEdgeID)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->setMeshTriangle(*meshID, *firstEdgeID, *secondEdgeID, *thirdEdgeID);
}

void precicef_set_triangle_we_(
    const int *meshID,
    const int *firstVertexID,
    const int *secondVertexID,
    const int *thirdVertexID)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->setMeshTriangleWithEdges(*meshID, *firstVertexID, *secondVertexID, *thirdVertexID);
}

void precicef_set_quad_(
    const int *meshID,
    const int *firstEdgeID,
    const int *secondEdgeID,
    const int *thirdEdgeID,
    const int *fourthEdgeID)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->setMeshQuad(*meshID, *firstEdgeID, *secondEdgeID, *thirdEdgeID, *fourthEdgeID);
}

void precicef_set_quad_we_(
    const int *meshID,
    const int *firstVertexID,
    const int *secondVertexID,
    const int *thirdVertexID,
    const int *fourthVertexID)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->setMeshQuadWithEdges(*meshID, *firstVertexID, *secondVertexID, *thirdVertexID, *fourthVertexID);
}

void precicef_write_bvdata_(
    const int *dataID,
    const int *size,
    int *      valueIndices,
    double *   values)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->writeBlockVectorData(*dataID, *size, valueIndices, values);
}

void precicef_write_vdata_(
    const int *   dataID,
    const int *   valueIndex,
    const double *dataValue)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->writeVectorData(*dataID, *valueIndex, dataValue);
}

void precicef_write_bsdata_(
    const int *dataID,
    const int *size,
    int *      valueIndices,
    double *   values)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->writeBlockScalarData(*dataID, *size, valueIndices, values);
}

void precicef_write_sdata_(
    const int *   dataID,
    const int *   valueIndex,
    const double *dataValue)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->writeScalarData(*dataID, *valueIndex, *dataValue);
}

void precicef_read_bvdata_(
    const int *dataID,
    const int *size,
    int *      valueIndices,
    double *   values)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->readBlockVectorData(*dataID, *size, valueIndices, values);
}

void precicef_read_vdata_(
    const int *dataID,
    const int *valueIndex,
    double *   dataValue)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->readVectorData(*dataID, *valueIndex, dataValue);
}

void precicef_read_bsdata_(
    const int *dataID,
    const int *size,
    int *      valueIndices,
    double *   values)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->readBlockScalarData(*dataID, *size, valueIndices, values);
}

void precicef_read_sdata_(
    const int *dataID,
    const int *valueIndex,
    double *   dataValue)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->readScalarData(*dataID, *valueIndex, *dataValue);
}

void precicef_map_write_data_from_(
    const int *meshID)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->mapWriteDataFrom(*meshID);
}

void precicef_map_read_data_to_(
    const int *meshID)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->mapReadDataTo(*meshID);
}

int precice::impl::strippedLength(
    const char *string,
    int         length)
{
  int i = length - 1;
  while (((string[i] == ' ') || (string[i] == 0)) && (i >= 0)) {
    i--;
  }
  return i + 1;
}

void precicef_action_write_iter_checkp_(
    char *nameAction,
    int   lengthNameAction)
{
  const std::string &name = precice::constants::actionWriteIterationCheckpoint();
  PRECICE_ASSERT(name.size() < (size_t) lengthNameAction, name.size(), lengthNameAction);
  for (size_t i = 0; i < name.size(); i++) {
    nameAction[i] = name[i];
  }
}

void precicef_action_write_initial_data_(
    char *nameAction,
    int   lengthNameAction)
{
  const std::string &name = precice::constants::actionWriteInitialData();
  PRECICE_ASSERT(name.size() < (size_t) lengthNameAction, name.size(), lengthNameAction);
  for (size_t i = 0; i < name.size(); i++) {
    nameAction[i] = name[i];
  }
}

void precicef_action_read_iter_checkp_(
    char *nameAction,
    int   lengthNameAction)
{
  const std::string &name = precice::constants::actionReadIterationCheckpoint();
  PRECICE_ASSERT(name.size() < (size_t) lengthNameAction, name.size(), lengthNameAction);
  for (size_t i = 0; i < name.size(); i++) {
    nameAction[i] = name[i];
  }
}

void precicef_get_version_information_(
    char *versionInfo,
    int   lengthVersionInfo)
{
  const std::string &versionInformation = precice::versionInformation;
  PRECICE_ASSERT(versionInformation.size() < (size_t) lengthVersionInfo, versionInformation.size(), lengthVersionInfo);
  for (size_t i = 0; i < versionInformation.size(); i++) {
    versionInfo[i] = versionInformation[i];
  }
}

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif
#ifdef __clang__
#pragma clang diagnostic pop
#endif
