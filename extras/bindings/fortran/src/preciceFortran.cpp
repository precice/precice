#include "precice/preciceFortran.hpp"
#include <cstddef>
#include <iostream>
#include <memory>
#ifndef PRECICE_NO_MPI
#include <mpi.h>
#endif
#include <stdexcept>
#include <string>
#include <string_view>
#include "logging/LogMacros.hpp"
#include "logging/Logger.hpp"
#include "precice/impl/versions.hpp"
#include "precice/precice.hpp"
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

static std::unique_ptr<precice::Participant> impl = nullptr;

static precice::logging::Logger _log("preciceFortran");

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
    const int  *solverProcessIndex,
    const int  *solverProcessSize,
    int         participantNameLength,
    int         configFileNameLength)
try {
  impl = std::make_unique<precice::Participant>(
      precice::impl::strippedStringView(participantName, participantNameLength),
      precice::impl::strippedStringView(configFileName, configFileNameLength),
      *solverProcessIndex,
      *solverProcessSize);
} catch (::precice::Error &e) {
  std::abort();
}

void precicef_create_with_communicator_(
    const char *participantName,
    const char *configFileName,
    const int  *solverProcessIndex,
    const int  *solverProcessSize,
    const int  *communicator,
    int         participantNameLength,
    int         configFileNameLength)
try {
#ifndef PRECICE_NO_MPI
  // MPI_Comm_f2c requires an MPI_Fint argument, which may differ from int.
  // However, the Fortran standard guarantees interoperability of the int type,
  // so that is what is passed. Do the conversion in case int is not identical
  // to MPI_Fint.
  MPI_Fint f_communicator = static_cast<MPI_Fint>(*communicator);
  auto     c_communicator = MPI_Comm_f2c(f_communicator);
  impl                    = std::make_unique<precice::Participant>(
      precice::impl::strippedStringView(participantName, participantNameLength),
      precice::impl::strippedStringView(configFileName, configFileNameLength),
      *solverProcessIndex,
      *solverProcessSize,
      &c_communicator);
#else
  PRECICE_WARN("preCICE was configured without MPI but you passed an MPI communicator. preCICE ignores the communicator and continues.");
  impl = std::make_unique<precice::Participant>(
      precice::impl::strippedStringView(participantName, participantNameLength),
      precice::impl::strippedStringView(configFileName, configFileNameLength),
      *solverProcessIndex,
      *solverProcessSize);
#endif
} catch (::precice::Error &e) {
  std::abort();
}

void precicef_initialize_()
try {
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->initialize();
} catch (::precice::Error &e) {
  std::abort();
}

void precicef_advance_(
    const double *timeStepSize)
try {
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->advance(*timeStepSize);
} catch (::precice::Error &e) {
  std::abort();
}

void precicef_finalize_()
try {
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->finalize();
  impl.reset();
} catch (::precice::Error &e) {
  std::abort();
}

void precicef_get_mesh_dimensions_(
    const char *meshName,
    int        *dimensions,
    int         meshNameLength)
try {
  PRECICE_CHECK(impl != nullptr, errormsg);
  *dimensions = impl->getMeshDimensions(precice::impl::strippedStringView(meshName, meshNameLength));
} catch (::precice::Error &e) {
  std::abort();
}

void precicef_get_data_dimensions_(
    const char *meshName,
    const char *dataName,
    int        *dimensions,
    int         meshNameLength,
    int         dataNameLength)
try {
  PRECICE_CHECK(impl != nullptr, errormsg);
  *dimensions = impl->getDataDimensions(precice::impl::strippedStringView(meshName, meshNameLength), precice::impl::strippedStringView(dataName, dataNameLength));
} catch (::precice::Error &e) {
  std::abort();
}

void precicef_is_coupling_ongoing_(
    int *isOngoing)
try {
  PRECICE_CHECK(impl != nullptr, errormsg);
  if (impl->isCouplingOngoing()) {
    *isOngoing = 1;
  } else {
    *isOngoing = 0;
  }
} catch (::precice::Error &e) {
  std::abort();
}

void precicef_is_time_window_complete_(
    int *isComplete)
try {
  PRECICE_CHECK(impl != nullptr, errormsg);
  if (impl->isTimeWindowComplete()) {
    *isComplete = 1;
  } else {
    *isComplete = 0;
  }
} catch (::precice::Error &e) {
  std::abort();
}

void precicef_get_max_time_step_size_(
    double *maxTimeStepSize)
try {
  PRECICE_CHECK(impl != nullptr, errormsg);
  *maxTimeStepSize = impl->getMaxTimeStepSize();
} catch (::precice::Error &e) {
  std::abort();
}

void precicef_requires_initial_data_(
    int *isRequired)
try {
  PRECICE_CHECK(impl != nullptr, errormsg);
  *isRequired = impl->requiresInitialData() ? 1 : 0;
} catch (::precice::Error &e) {
  std::abort();
}

void precicef_requires_writing_checkpoint_(
    int *isRequired)
try {
  PRECICE_CHECK(impl != nullptr, errormsg);
  *isRequired = impl->requiresWritingCheckpoint() ? 1 : 0;
} catch (::precice::Error &e) {
  std::abort();
}

void precicef_requires_reading_checkpoint_(
    int *isRequired)
try {
  PRECICE_CHECK(impl != nullptr, errormsg);
  *isRequired = impl->requiresReadingCheckpoint() ? 1 : 0;
} catch (::precice::Error &e) {
  std::abort();
}

void precicef_requires_mesh_connectivity_for_(
    const char *meshName,
    int        *required,
    int         meshNameLength)
try {
  PRECICE_CHECK(impl != nullptr, errormsg);
  if (impl->requiresMeshConnectivityFor(precice::impl::strippedStringView(meshName, meshNameLength))) {
    *required = 1;
  } else {
    *required = 0;
  }
} catch (::precice::Error &e) {
  std::abort();
}

void precicef_reset_mesh_(
    const char *meshName,
    int         meshNameLength)
try {
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->resetMesh(precice::impl::strippedStringView(meshName, meshNameLength));
} catch (::precice::Error &e) {
  std::abort();
}

void precicef_set_vertex_(
    const char   *meshName,
    const double *position,
    int          *vertexID,
    int           meshNameLength)
try {
  PRECICE_CHECK(impl != nullptr, errormsg);
  auto sv           = precice::impl::strippedStringView(meshName, meshNameLength);
  auto positionSize = static_cast<unsigned long>(impl->getMeshDimensions(sv));
  *vertexID         = impl->setMeshVertex(sv, {position, positionSize});
} catch (::precice::Error &e) {
  std::abort();
}

void precicef_get_mesh_vertex_size_(
    const char *meshName,
    int        *meshSize,
    int         meshNameLength)
try {
  PRECICE_CHECK(impl != nullptr, errormsg);
  *meshSize = impl->getMeshVertexSize(precice::impl::strippedStringView(meshName, meshNameLength));
} catch (::precice::Error &e) {
  std::abort();
}

void precicef_set_vertices_(
    const char *meshName,
    const int  *size,
    double     *coordinates,
    int        *ids,
    int         meshNameLength)
try {
  PRECICE_CHECK(impl != nullptr, errormsg);
  auto sv           = precice::impl::strippedStringView(meshName, meshNameLength);
  auto positionSize = static_cast<unsigned long>(impl->getMeshDimensions(sv) * *size);
  impl->setMeshVertices(sv, {coordinates, positionSize}, {ids, static_cast<unsigned long>(*size)});
} catch (::precice::Error &e) {
  std::abort();
}

void precicef_set_edge_(
    const char *meshName,
    const int  *firstVertexID,
    const int  *secondVertexID,
    int         meshNameLength)
try {
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->setMeshEdge(precice::impl::strippedStringView(meshName, meshNameLength), *firstVertexID, *secondVertexID);
} catch (::precice::Error &e) {
  std::abort();
}

void precicef_set_mesh_edges_(
    const char *meshName,
    const int  *size,
    const int  *ids,
    int         meshNameLength)
try {
  PRECICE_CHECK(impl != nullptr, errormsg);
  auto idsSize = static_cast<unsigned long>(*size) * 2;
  impl->setMeshEdges(precice::impl::strippedStringView(meshName, meshNameLength), {ids, idsSize});
} catch (::precice::Error &e) {
  std::abort();
}

void precicef_set_triangle_(
    const char *meshName,
    const int  *firstVertexID,
    const int  *secondVertexID,
    const int  *thirdVertexID,
    int         meshNameLength)
try {
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->setMeshTriangle(precice::impl::strippedStringView(meshName, meshNameLength), *firstVertexID, *secondVertexID, *thirdVertexID);
} catch (::precice::Error &e) {
  std::abort();
}

void precicef_set_mesh_triangles_(
    const char *meshName,
    const int  *size,
    const int  *ids,
    int         meshNameLength)
try {
  PRECICE_CHECK(impl != nullptr, errormsg);
  auto idsSize = static_cast<unsigned long>(*size) * 3;
  impl->setMeshTriangles(precice::impl::strippedStringView(meshName, meshNameLength), {ids, idsSize});
} catch (::precice::Error &e) {
  std::abort();
}

void precicef_set_quad_(
    const char *meshName,
    const int  *firstVertexID,
    const int  *secondVertexID,
    const int  *thirdVertexID,
    const int  *fourthVertexID,
    int         meshNameLength)
try {
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->setMeshQuad(precice::impl::strippedStringView(meshName, meshNameLength), *firstVertexID, *secondVertexID, *thirdVertexID, *fourthVertexID);
} catch (::precice::Error &e) {
  std::abort();
}

void precicef_set_mesh_quads_(
    const char *meshName,
    const int  *size,
    const int  *ids,
    int         meshNameLength)
try {
  PRECICE_CHECK(impl != nullptr, errormsg);
  auto idsSize = static_cast<unsigned long>(*size) * 4;
  impl->setMeshQuads(precice::impl::strippedStringView(meshName, meshNameLength), {ids, idsSize});
} catch (::precice::Error &e) {
  std::abort();
}

void precicef_set_tetrahedron(
    const char *meshName,
    const int  *firstVertexID,
    const int  *secondVertexID,
    const int  *thirdVertexID,
    const int  *fourthVertexID,
    int         meshNameLength)
try {
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->setMeshTetrahedron(precice::impl::strippedStringView(meshName, meshNameLength), *firstVertexID, *secondVertexID, *thirdVertexID, *fourthVertexID);
} catch (::precice::Error &e) {
  std::abort();
}

void precicef_set_mesh_tetrahedra_(
    const char *meshName,
    const int  *size,
    const int  *ids,
    int         meshNameLength)
try {
  PRECICE_CHECK(impl != nullptr, errormsg);
  auto idsSize = static_cast<unsigned long>(*size) * 4;
  impl->setMeshTetrahedra(precice::impl::strippedStringView(meshName, meshNameLength), {ids, idsSize});
} catch (::precice::Error &e) {
  std::abort();
}

void precicef_write_data_(
    const char *meshName,
    const char *dataName,
    const int  *size,
    int        *ids,
    double     *values,
    int         meshNameLength,
    int         dataNameLength)
try {
  PRECICE_CHECK(impl != nullptr, errormsg);
  auto strippedMeshName = precice::impl::strippedStringView(meshName, meshNameLength);
  auto strippedDataName = precice::impl::strippedStringView(dataName, dataNameLength);
  auto dataSize         = *size * impl->getDataDimensions(strippedMeshName, strippedDataName);
  impl->writeData(strippedMeshName,
                  strippedDataName,
                  {ids, static_cast<unsigned long>(*size)},
                  {values, static_cast<unsigned long>(dataSize)});
} catch (::precice::Error &e) {
  std::abort();
}

void precicef_read_data_(
    const char   *meshName,
    const char   *dataName,
    const int    *size,
    int          *ids,
    const double *relativeReadTime,
    double       *values,
    int           meshNameLength,
    int           dataNameLength)
try {
  PRECICE_CHECK(impl != nullptr, errormsg);
  auto strippedMeshName = precice::impl::strippedStringView(meshName, meshNameLength);
  auto strippedDataName = precice::impl::strippedStringView(dataName, dataNameLength);
  auto dataSize         = *size * impl->getDataDimensions(strippedMeshName, strippedDataName);
  impl->readData(
      strippedMeshName,
      strippedDataName,
      {ids, static_cast<unsigned long>(*size)},
      *relativeReadTime,
      {values, static_cast<unsigned long>(dataSize)});
} catch (::precice::Error &e) {
  std::abort();
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
  std::strncpy(versionInfo, precice::versionInformation, lengthVersionInfo);
}

void precicef_requires_gradient_data_for_(
    const char *meshName,
    const char *dataName, int *required,
    int meshNameLength,
    int dataNameLength)
try {
  PRECICE_CHECK(impl != nullptr, errormsg);
  if (impl->requiresGradientDataFor(precice::impl::strippedStringView(meshName, meshNameLength), precice::impl::strippedStringView(dataName, dataNameLength))) {
    *required = 1;
  } else {
    *required = 0;
  }
} catch (::precice::Error &e) {
  std::abort();
}

void precicef_write_gradient_data_(
    const char   *meshName,
    const char   *dataName,
    const int    *size,
    const int    *ids,
    const double *gradients,
    int           meshNameLength,
    int           dataNameLength)
try {
  PRECICE_CHECK(impl != nullptr, errormsg);
  auto strippedMeshName   = precice::impl::strippedStringView(meshName, meshNameLength);
  auto strippedDataName   = precice::impl::strippedStringView(dataName, dataNameLength);
  auto gradientComponents = impl->getMeshDimensions(strippedMeshName) * impl->getDataDimensions(strippedMeshName, strippedDataName);
  auto gradientSize       = *size * gradientComponents;

  impl->writeGradientData(strippedMeshName,
                          strippedDataName,
                          {ids, static_cast<unsigned long>(*size)},
                          {gradients, static_cast<unsigned long>(gradientSize)});
} catch (::precice::Error &e) {
  std::abort();
}

void precicef_write_and_map_data_(
    const char *meshName,
    const char *dataName,
    const int  *size,
    double     *coordinates,
    double     *values,
    int         meshNameLength,
    int         dataNameLength)
try {
  PRECICE_CHECK(impl != nullptr, errormsg);
  auto strippedMeshName = precice::impl::strippedStringView(meshName, meshNameLength);
  auto strippedDataName = precice::impl::strippedStringView(dataName, dataNameLength);
  auto coordinatesSize  = *size * impl->getMeshDimensions(strippedMeshName);
  auto dataSize         = *size * impl->getDataDimensions(strippedMeshName, strippedDataName);
  impl->writeAndMapData(strippedMeshName,
                        strippedDataName,
                        {coordinates, static_cast<unsigned long>(coordinatesSize)},
                        {values, static_cast<unsigned long>(dataSize)});
} catch (::precice::Error &e) {
  std::abort();
}

void precicef_map_and_read_data_(
    const char   *meshName,
    const char   *dataName,
    const int    *size,
    double       *coordinates,
    const double *relativeReadTime,
    double       *values,
    int           meshNameLength,
    int           dataNameLength)
try {
  PRECICE_CHECK(impl != nullptr, errormsg);
  auto strippedMeshName = precice::impl::strippedStringView(meshName, meshNameLength);
  auto strippedDataName = precice::impl::strippedStringView(dataName, dataNameLength);
  auto coordinatesSize  = *size * impl->getMeshDimensions(strippedMeshName);
  auto dataSize         = *size * impl->getDataDimensions(strippedMeshName, strippedDataName);
  impl->mapAndReadData(
      strippedMeshName,
      strippedDataName,
      {coordinates, static_cast<unsigned long>(coordinatesSize)},
      *relativeReadTime,
      {values, static_cast<unsigned long>(dataSize)});
} catch (::precice::Error &e) {
  std::abort();
}

void precicef_set_mesh_access_region_(
    const char   *meshName,
    const double *boundingBox,
    int           meshNameLength)
try {
  PRECICE_CHECK(impl != nullptr, errormsg);
  auto sv     = precice::impl::strippedStringView(meshName, meshNameLength);
  auto bbSize = static_cast<unsigned long>(impl->getMeshDimensions(sv) * 2);
  impl->setMeshAccessRegion(sv, {boundingBox, bbSize});
} catch (::precice::Error &e) {
  std::abort();
}

void precicef_get_mesh_vertex_ids_and_coordinates_(
    const char *meshName,
    const int  *size,
    int        *ids,
    double     *coordinates,
    int         meshNameLength)
try {
  PRECICE_CHECK(impl != nullptr, errormsg);
  auto sv              = precice::impl::strippedStringView(meshName, meshNameLength);
  auto coordinatesSize = static_cast<unsigned long>(impl->getMeshDimensions(sv) * *size);
  impl->getMeshVertexIDsAndCoordinates(sv, {ids, static_cast<unsigned long>(*size)}, {coordinates, coordinatesSize});
} catch (::precice::Error &e) {
  std::abort();
}

void precicef_start_profiling_section_(
    const char *sectionName,
    int         sectionNameLength)
try {
  PRECICE_CHECK(impl != nullptr, errormsg);
  auto sv = precice::impl::strippedStringView(sectionName, sectionNameLength);
  impl->startProfilingSection(sv);
} catch (::precice::Error &e) {
  std::abort();
}

void precicef_stop_last_profiling_section_()
try {
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->stopLastProfilingSection();
} catch (::precice::Error &e) {
  std::abort();
}

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif
#ifdef __clang__
#pragma clang diagnostic pop
#endif
