#include "precice/SolverInterfaceFortran.hpp"
#include <iostream>
#include <memory>
#include <stddef.h>
#include <stdexcept>
#include <string>
#include <string_view>
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

namespace precice::impl {
/**
 * @brief Returns length of string without trailing whitespace.
 */
int strippedLength(const char *string, int length);

std::string_view strippedStringView(const char *string, int length);

} // namespace precice::impl

void precicef_create_(
    const char *participantName,
    const char *configFileName,
    const int * solverProcessIndex,
    const int * solverProcessSize,
    int         participantNameLength,
    int         configFileNameLength)
{
  impl.reset(new precice::SolverInterface(
      precice::impl::strippedStringView(participantName, participantNameLength),
      precice::impl::strippedStringView(configFileName, configFileNameLength),
      *solverProcessIndex,
      *solverProcessSize));
}

void precicef_initialize_()
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->initialize();
}

void precicef_advance_(
    const double *timeStepSize)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->advance(*timeStepSize);
}

void precicef_finalize_()
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->finalize();
  impl.reset();
}

void precicef_get_mesh_dimensions_(
    const char *meshName,
    int *       dimensions,
    int         meshNameLength)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  *dimensions = impl->getMeshDimensions(precice::impl::strippedStringView(meshName, meshNameLength));
}

void precicef_get_data_dimensions_(
    const char *meshName,
    const char *dataName,
    int *       dimensions,
    int         meshNameLength,
    int         dataNameLength)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  *dimensions = impl->getDataDimensions(precice::impl::strippedStringView(meshName, meshNameLength), precice::impl::strippedStringView(dataName, dataNameLength));
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

void precicef_get_max_time_step_size_(
    double *maxTimeStepSize)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  *maxTimeStepSize = impl->getMaxTimeStepSize();
}

void precicef_requires_initial_data_(
    int *isRequired)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  *isRequired = impl->requiresInitialData() ? 1 : 0;
}

void precicef_requires_writing_checkpoint_(
    int *isRequired)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  *isRequired = impl->requiresWritingCheckpoint() ? 1 : 0;
}

void precicef_requires_reading_checkpoint_(
    int *isRequired)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  *isRequired = impl->requiresReadingCheckpoint() ? 1 : 0;
}

void precicef_has_mesh_(
    const char *meshName,
    int *       hasMesh,
    int         meshNameLength)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  if (impl->hasMesh(precice::impl::strippedStringView(meshName, meshNameLength))) {
    *hasMesh = 1;
  } else {
    *hasMesh = 0;
  }
}

void precicef_has_data_(
    const char *meshName,
    const char *dataName,
    int *       hasData,
    int         meshNameLength,
    int         dataNameLength)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  if (impl->hasData(precice::impl::strippedStringView(meshName, meshNameLength), precice::impl::strippedStringView(dataName, dataNameLength))) {
    *hasData = 1;
  } else {
    *hasData = 0;
  }
}

void precicef_requires_mesh_connectivity_for_(
    const char *meshName,
    int *       required,
    int         meshNameLength)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  if (impl->requiresMeshConnectivityFor(precice::impl::strippedStringView(meshName, meshNameLength))) {
    *required = 1;
  } else {
    *required = 0;
  }
}

void precicef_set_vertex_(
    const char *  meshName,
    const double *position,
    int *         vertexID,
    int           meshNameLength)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  *vertexID = impl->setMeshVertex(precice::impl::strippedStringView(meshName, meshNameLength), position);
}

void precicef_get_mesh_vertex_size_(
    const char *meshName,
    int *       meshSize,
    int         meshNameLength)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  *meshSize = impl->getMeshVertexSize(precice::impl::strippedStringView(meshName, meshNameLength));
}

void precicef_set_vertices_(
    const char *meshName,
    const int * size,
    double *    positions,
    int *       positionIDs,
    int         meshNameLength)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->setMeshVertices(precice::impl::strippedStringView(meshName, meshNameLength), *size, positions, positionIDs);
}

void precicef_set_edge_(
    const char *meshName,
    const int * firstVertexID,
    const int * secondVertexID,
    int         meshNameLength)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->setMeshEdge(precice::impl::strippedStringView(meshName, meshNameLength), *firstVertexID, *secondVertexID);
}

void precicef_set_mesh_edges_(
    const char *meshName,
    const int * size,
    const int * vertices,
    int         meshNameLength)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->setMeshEdges(precice::impl::strippedStringView(meshName, meshNameLength), *size, vertices);
}

void precicef_set_triangle_(
    const char *meshName,
    const int * firstVertexID,
    const int * secondVertexID,
    const int * thirdVertexID,
    int         meshNameLength)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->setMeshTriangle(precice::impl::strippedStringView(meshName, meshNameLength), *firstVertexID, *secondVertexID, *thirdVertexID);
}

void precicef_set_mesh_triangles_(
    const char *meshName,
    const int * size,
    const int * vertices,
    int         meshNameLength)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->setMeshTriangles(precice::impl::strippedStringView(meshName, meshNameLength), *size, vertices);
}

void precicef_set_quad_(
    const char *meshName,
    const int * firstVertexID,
    const int * secondVertexID,
    const int * thirdVertexID,
    const int * fourthVertexID,
    int         meshNameLength)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->setMeshQuad(precice::impl::strippedStringView(meshName, meshNameLength), *firstVertexID, *secondVertexID, *thirdVertexID, *fourthVertexID);
}

void precicef_set_mesh_quads_(
    const char *meshName,
    const int * size,
    const int * vertices,
    int         meshNameLength)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->setMeshQuads(precice::impl::strippedStringView(meshName, meshNameLength), *size, vertices);
}

void precicef_set_tetrahedron(
    const char *meshName,
    const int * firstVertexID,
    const int * secondVertexID,
    const int * thirdVertexID,
    const int * fourthVertexID,
    int         meshNameLength)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->setMeshTetrahedron(precice::impl::strippedStringView(meshName, meshNameLength), *firstVertexID, *secondVertexID, *thirdVertexID, *fourthVertexID);
}

void precicef_set_mesh_tetrahedra_(
    const char *meshName,
    const int * size,
    const int * vertices,
    int         meshNameLength)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->setMeshTetrahedra(precice::impl::strippedStringView(meshName, meshNameLength), *size, vertices);
}

void precicef_write_bvdata_(
    const char *meshName,
    const char *dataName,
    const int * size,
    int *       valueIndices,
    double *    values,
    int         meshNameLength,
    int         dataNameLength)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->writeBlockVectorData(precice::impl::strippedStringView(meshName, meshNameLength), precice::impl::strippedStringView(dataName, dataNameLength), *size, valueIndices, values);
}

void precicef_write_vdata_(
    const char *  meshName,
    const char *  dataName,
    const int *   valueIndex,
    const double *dataValue,
    int           meshNameLength,
    int           dataNameLength)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->writeVectorData(precice::impl::strippedStringView(meshName, meshNameLength), precice::impl::strippedStringView(dataName, dataNameLength), *valueIndex, dataValue);
}

void precicef_write_bsdata_(
    const char *meshName,
    const char *dataName,
    const int * size,
    int *       valueIndices,
    double *    values,
    int         meshNameLength,
    int         dataNameLength)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->writeBlockScalarData(precice::impl::strippedStringView(meshName, meshNameLength), precice::impl::strippedStringView(dataName, dataNameLength), *size, valueIndices, values);
}

void precicef_write_sdata_(
    const char *  meshName,
    const char *  dataName,
    const int *   valueIndex,
    const double *dataValue,
    int           meshNameLength,
    int           dataNameLength)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->writeScalarData(precice::impl::strippedStringView(meshName, meshNameLength), precice::impl::strippedStringView(dataName, dataNameLength), *valueIndex, *dataValue);
}

void precicef_read_bvdata_(
    const char *  meshName,
    const char *  dataName,
    const int *   size,
    int *         valueIndices,
    const double *relativeReadTime,
    double *      values,
    int           meshNameLength,
    int           dataNameLength)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->readBlockVectorData(precice::impl::strippedStringView(meshName, meshNameLength), precice::impl::strippedStringView(dataName, dataNameLength), *size, valueIndices, *relativeReadTime, values);
}

void precicef_read_vdata_(
    const char *  meshName,
    const char *  dataName,
    const int *   valueIndex,
    const double *relativeReadTime,
    double *      dataValue,
    int           meshNameLength,
    int           dataNameLength)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->readVectorData(precice::impl::strippedStringView(meshName, meshNameLength), precice::impl::strippedStringView(dataName, dataNameLength), *valueIndex, *relativeReadTime, dataValue);
}

void precicef_read_bsdata_(
    const char *  meshName,
    const char *  dataName,
    const int *   size,
    int *         valueIndices,
    const double *relativeReadTime,
    double *      values,
    int           meshNameLength,
    int           dataNameLength)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->readBlockScalarData(precice::impl::strippedStringView(meshName, meshNameLength), precice::impl::strippedStringView(dataName, dataNameLength), *size, valueIndices, *relativeReadTime, values);
}

void precicef_read_sdata_(
    const char *  meshName,
    const char *  dataName,
    const int *   valueIndex,
    const double *relativeReadTime,
    double *      dataValue,
    int           meshNameLength,
    int           dataNameLength)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->readScalarData(precice::impl::strippedStringView(meshName, meshNameLength), precice::impl::strippedStringView(dataName, dataNameLength), *valueIndex, *relativeReadTime, *dataValue);
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

std::string_view precice::impl::strippedStringView(const char *string, int length)
{
  return {string, static_cast<std::string_view::size_type>(strippedLength(string, length))};
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

void precicef_requires_gradient_data_for_(
    const char *meshName,
    const char *dataName, int *required,
    int meshNameLength,
    int dataNameLength)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  if (impl->requiresGradientDataFor(precice::impl::strippedStringView(meshName, meshNameLength), precice::impl::strippedStringView(dataName, dataNameLength))) {
    *required = 1;
  } else {
    *required = 0;
  }
}

void precicef_write_sgradient_data_(
    const char *  meshName,
    const char *  dataName,
    const int *   valueIndex,
    const double *gradientValues,
    int           meshNameLength,
    int           dataNameLength)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->writeScalarGradientData(precice::impl::strippedStringView(meshName, meshNameLength), precice::impl::strippedStringView(dataName, dataNameLength), *valueIndex, gradientValues);
}

void precicef_write_bsgradient_data_(
    const char *  meshName,
    const char *  dataName,
    const int *   size,
    const int *   valueIndices,
    const double *gradientValues,
    int           meshNameLength,
    int           dataNameLength)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->writeBlockScalarGradientData(precice::impl::strippedStringView(meshName, meshNameLength), precice::impl::strippedStringView(dataName, dataNameLength), *size, valueIndices, gradientValues);
}

void precicef_write_vgradient_data_(
    const char *  meshName,
    const char *  dataName,
    const int *   valueIndex,
    const double *gradientValues,
    int           meshNameLength,
    int           dataNameLength)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->writeVectorGradientData(precice::impl::strippedStringView(meshName, meshNameLength), precice::impl::strippedStringView(dataName, dataNameLength), *valueIndex, gradientValues);
}

void precicef_write_bvgradient_data_(
    const char *  meshName,
    const char *  dataName,
    const int *   size,
    const int *   valueIndices,
    const double *gradientValues,
    int           meshNameLength,
    int           dataNameLength)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->writeBlockVectorGradientData(precice::impl::strippedStringView(meshName, meshNameLength), precice::impl::strippedStringView(dataName, dataNameLength), *size, valueIndices, gradientValues);
}

void precicef_set_mesh_access_region_(
    const char *  meshName,
    const double *boundingBox,
    int           meshNameLength)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->setMeshAccessRegion(precice::impl::strippedStringView(meshName, meshNameLength), boundingBox);
}

void precicef_get_mesh_vertices_and_ids_(
    const char *meshName,
    const int   size,
    int *       ids,
    double *    coordinates,
    int         meshNameLength)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->getMeshVerticesAndIDs(precice::impl::strippedStringView(meshName, meshNameLength), size, ids, coordinates);
}

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif
#ifdef __clang__
#pragma clang diagnostic pop
#endif
