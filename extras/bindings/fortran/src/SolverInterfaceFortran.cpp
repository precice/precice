#include "precice/SolverInterfaceFortran.hpp"
#include <iostream>
#include <memory>
#include <stddef.h>
#include <string>
#include "logging/LogMacros.hpp"
#include "logging/Logger.hpp"
#include "precice/SolverInterface.hpp"
#include "precice/exceptions.hpp"
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
using Err                   = ::precice::APIError;

namespace precice::impl {
/**
 * @brief Returns length of string without trailing whitespace.
 */
int strippedLength(const char *string, int length);
} // namespace precice::impl

void precicef_create_(
    const char *participantName,
    const char *configFileName,
    const int * solverProcessIndex,
    const int * solverProcessSize,
    int         lengthAccessorName,
    int         lengthConfigFileName)
{
  // cout << "lengthAccessorName: " << lengthAccessorName << '\n';
  // cout << "lengthConfigFileName: " << lengthConfigFileName << '\n';
  // cout << "solverProcessIndex: " << *solverProcessIndex << '\n';
  // cout << "solverProcessSize: " << *solverProcessSize << '\n';
  int    strippedLength = precice::impl::strippedLength(participantName, lengthAccessorName);
  string stringAccessorName(participantName, strippedLength);
  strippedLength = precice::impl::strippedLength(configFileName, lengthConfigFileName);
  string stringConfigFileName(configFileName, strippedLength);
  // cout << "Accessor: " << stringAccessorName << "!" << '\n';
  // cout << "Config  : " << stringConfigFileName << "!" << '\n';
  impl.reset(new precice::SolverInterface(stringAccessorName,
                                          stringConfigFileName,
                                          *solverProcessIndex, *solverProcessSize));
}

void precicef_initialize_(
    double *timestepLengthLimit)
{
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  *timestepLengthLimit = impl->initialize();
}

void precicef_advance_(
    double *timestepLengthLimit)
{
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  *timestepLengthLimit = impl->advance(*timestepLengthLimit);
}

void precicef_finalize_()
{
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  impl->finalize();
  impl.reset();
}

void precicef_get_dims_(
    int *dimensions)
{
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  *dimensions = impl->getDimensions();
}

void precicef_is_coupling_ongoing_(
    int *isOngoing)
{
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  if (impl->isCouplingOngoing()) {
    *isOngoing = 1;
  } else {
    *isOngoing = 0;
  }
}

void precicef_is_time_window_complete_(
    int *isComplete)
{
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  if (impl->isTimeWindowComplete()) {
    *isComplete = 1;
  } else {
    *isComplete = 0;
  }
}

void precicef_requires_initial_data_(
    int *isRequired)
{
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  *isRequired = impl->requiresInitialData() ? 1 : 0;
}

void precicef_requires_writing_checkpoint_(
    int *isRequired)
{
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  *isRequired = impl->requiresWritingCheckpoint() ? 1 : 0;
}

void precicef_requires_reading_checkpoint_(
    int *isRequired)
{
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  *isRequired = impl->requiresReadingCheckpoint() ? 1 : 0;
}

void precicef_has_mesh_(
    const char *meshName,
    int *       hasMesh,
    int         lengthMeshName)
{
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
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
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
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
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
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
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  int    strippedLength = precice::impl::strippedLength(dataName, lengthDataName);
  string stringDataName(dataName, strippedLength);
  *dataID = impl->getDataID(stringDataName, *meshID);
}

void precicef_requires_mesh_connectivity_for_(
    const int *meshID,
    int *      required)
{
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  if (impl->requiresMeshConnectivityFor(*meshID)) {
    *required = 1;
  } else {
    *required = 0;
  }
}

void precicef_set_vertex_(
    const int *   meshID,
    const double *position,
    int *         vertexID)
{
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  *vertexID = impl->setMeshVertex(*meshID, position);
}

void precicef_get_mesh_vertex_size_(
    const int *meshID,
    int *      meshSize)
{
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  *meshSize = impl->getMeshVertexSize(*meshID);
}

void precicef_set_vertices_(
    const int *meshID,
    const int *size,
    double *   positions,
    int *      positionIDs)
{
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  impl->setMeshVertices(*meshID, *size, positions, positionIDs);
}

void precicef_set_edge_(
    const int *meshID,
    const int *firstVertexID,
    const int *secondVertexID)
{
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  impl->setMeshEdge(*meshID, *firstVertexID, *secondVertexID);
}

void precicef_set_mesh_edges_(
    const int *meshID,
    const int *size,
    const int *vertices)
{
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  impl->setMeshEdges(*meshID, *size, vertices);
}

void precicef_set_triangle_(
    const int *meshID,
    const int *firstVertexID,
    const int *secondVertexID,
    const int *thirdVertexID)
{
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  impl->setMeshTriangle(*meshID, *firstVertexID, *secondVertexID, *thirdVertexID);
}

void precicef_set_mesh_triangles_(
    const int *meshID,
    const int *size,
    const int *vertices)
{
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  impl->setMeshTriangles(*meshID, *size, vertices);
}

void precicef_set_quad_(
    const int *meshID,
    const int *firstVertexID,
    const int *secondVertexID,
    const int *thirdVertexID,
    const int *fourthVertexID)
{
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  impl->setMeshQuad(*meshID, *firstVertexID, *secondVertexID, *thirdVertexID, *fourthVertexID);
}

void precicef_set_mesh_quads_(
    const int *meshID,
    const int *size,
    const int *vertices)
{
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  impl->setMeshQuads(*meshID, *size, vertices);
}

void precicef_set_tetrahedron(
    const int *meshID,
    const int *firstVertexID,
    const int *secondVertexID,
    const int *thirdVertexID,
    const int *fourthVertexID)
{
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  impl->setMeshTetrahedron(*meshID, *firstVertexID, *secondVertexID, *thirdVertexID, *fourthVertexID);
}

void precicef_set_mesh_tetrahedra_(
    const int *meshID,
    const int *size,
    const int *vertices)
{
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  impl->setMeshTetrahedra(*meshID, *size, vertices);
}

void precicef_write_bvdata_(
    const int *dataID,
    const int *size,
    int *      valueIndices,
    double *   values)
{
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  impl->writeBlockVectorData(*dataID, *size, valueIndices, values);
}

void precicef_write_vdata_(
    const int *   dataID,
    const int *   valueIndex,
    const double *dataValue)
{
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  impl->writeVectorData(*dataID, *valueIndex, dataValue);
}

void precicef_write_bsdata_(
    const int *dataID,
    const int *size,
    int *      valueIndices,
    double *   values)
{
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  impl->writeBlockScalarData(*dataID, *size, valueIndices, values);
}

void precicef_write_sdata_(
    const int *   dataID,
    const int *   valueIndex,
    const double *dataValue)
{
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  impl->writeScalarData(*dataID, *valueIndex, *dataValue);
}

void precicef_read_bvdata_(
    const int *dataID,
    const int *size,
    int *      valueIndices,
    double *   values)
{
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  impl->readBlockVectorData(*dataID, *size, valueIndices, values);
}

void precicef_read_vdata_(
    const int *dataID,
    const int *valueIndex,
    double *   dataValue)
{
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  impl->readVectorData(*dataID, *valueIndex, dataValue);
}

void precicef_read_bsdata_(
    const int *dataID,
    const int *size,
    int *      valueIndices,
    double *   values)
{
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  impl->readBlockScalarData(*dataID, *size, valueIndices, values);
}

void precicef_read_sdata_(
    const int *dataID,
    const int *valueIndex,
    double *   dataValue)
{
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  impl->readScalarData(*dataID, *valueIndex, *dataValue);
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

void precicef_set_mesh_access_region_(
    const int     meshID,
    const double *boundingBox)
{
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  impl->setMeshAccessRegion(meshID, boundingBox);
}

void precicef_get_mesh_vertices_and_IDs_(
    const int meshID,
    const int size,
    int *     ids,
    double *  coordinates)
{
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  impl->getMeshVerticesAndIDs(meshID, size, ids, coordinates);
}

void precicef_requires_gradient_data_for_(const int *dataID, int *required)
{
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  if (impl->requiresGradientDataFor(*dataID)) {
    *required = 1;
  } else {
    *required = 0;
  }
}

void precicef_write_sgradient_data_(
    const int *   dataID,
    const int *   valueIndex,
    const double *gradientValues)
{
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  impl->writeScalarGradientData(*dataID, *valueIndex, gradientValues);
}

void precicef_write_bsgradient_data_(
    const int *   dataID,
    const int *   size,
    const int *   valueIndices,
    const double *gradientValues)
{
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  impl->writeBlockScalarGradientData(*dataID, *size, valueIndices, gradientValues);
}

void precicef_write_vgradient_data_(
    const int *   dataID,
    const int *   valueIndex,
    const double *gradientValues)
{
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  impl->writeVectorGradientData(*dataID, *valueIndex, gradientValues);
}

void precicef_write_bvgradient_data_(
    const int *   dataID,
    const int *   size,
    const int *   valueIndices,
    const double *gradientValues)
{
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  impl->writeBlockVectorGradientData(*dataID, *size, valueIndices, gradientValues);
}

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif
#ifdef __clang__
#pragma clang diagnostic pop
#endif
