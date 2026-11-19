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
#include "precice/impl/versions.hpp"
#include "precice/precice.hpp"

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
#endif

using namespace std;

namespace precice::bindings::fortran {

static std::unique_ptr<precice::Participant> impl = nullptr;

/// Returns length of string without trailing whitespace.
int strippedLength(const char *string, int length)
{
  int i = length - 1;
  while (((string[i] == ' ') || (string[i] == 0)) && (i >= 0)) {
    i--;
  }
  return i + 1;
}

std::string_view strippedStringView(const char *string, int length)
{
  return {string, static_cast<std::string_view::size_type>(strippedLength(string, length))};
}

void ensureNotConstructed()
{
  if (impl != nullptr) {
    std::cerr << "ERROR: The participant has already been created. Be sure to call \"precicef_finalize\" first.\n";
    std::abort();
  }
}

void ensureConstructed()
{
  if (impl == nullptr) {
    std::cerr << "ERROR: The participant has not yet been created. Be sure to call \"precicef_create\" first.\n";
    std::abort();
  }
}

} // namespace precice::bindings::fortran

using namespace precice::bindings::fortran;

void precicef_create_(
    const char *participantName,
    const char *configFileName,
    const int  *solverProcessIndex,
    const int  *solverProcessSize,
    int         participantNameLength,
    int         configFileNameLength)
try {
  ensureNotConstructed();
  impl = std::make_unique<precice::Participant>(
      strippedStringView(participantName, participantNameLength),
      strippedStringView(configFileName, configFileNameLength),
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
  ensureNotConstructed();
#ifndef PRECICE_NO_MPI
  // MPI_Comm_f2c requires an MPI_Fint argument, which may differ from int.
  // However, the Fortran standard guarantees interoperability of the int type,
  // so that is what is passed. Do the conversion in case int is not identical
  // to MPI_Fint.
  MPI_Fint f_communicator = static_cast<MPI_Fint>(*communicator);
  auto     c_communicator = MPI_Comm_f2c(f_communicator);
  impl                    = std::make_unique<precice::Participant>(
      strippedStringView(participantName, participantNameLength),
      strippedStringView(configFileName, configFileNameLength),
      *solverProcessIndex,
      *solverProcessSize,
      &c_communicator);
#else
  std::cerr << "WARNING: preCICE was configured without MPI but you passed an MPI communicator. preCICE ignores the communicator and continues.";
  impl = std::make_unique<precice::Participant>(
      strippedStringView(participantName, participantNameLength),
      strippedStringView(configFileName, configFileNameLength),
      *solverProcessIndex,
      *solverProcessSize);
#endif
} catch (::precice::Error &e) {
  std::abort();
}

void precicef_initialize_()
try {
  ensureConstructed();
  impl->initialize();
} catch (::precice::Error &e) {
  std::abort();
}

void precicef_advance_(
    const double *timeStepSize)
try {
  ensureConstructed();
  impl->advance(*timeStepSize);
} catch (::precice::Error &e) {
  std::abort();
}

void precicef_finalize_()
try {
  ensureConstructed();
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
  ensureConstructed();
  *dimensions = impl->getMeshDimensions(strippedStringView(meshName, meshNameLength));
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
  ensureConstructed();
  *dimensions = impl->getDataDimensions(strippedStringView(meshName, meshNameLength), strippedStringView(dataName, dataNameLength));
} catch (::precice::Error &e) {
  std::abort();
}

void precicef_is_coupling_ongoing_(
    int *isOngoing)
try {
  ensureConstructed();
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
  ensureConstructed();
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
  ensureConstructed();
  *maxTimeStepSize = impl->getMaxTimeStepSize();
} catch (::precice::Error &e) {
  std::abort();
}

void precicef_requires_initial_data_(
    int *isRequired)
try {
  ensureConstructed();
  *isRequired = impl->requiresInitialData() ? 1 : 0;
} catch (::precice::Error &e) {
  std::abort();
}

void precicef_requires_writing_checkpoint_(
    int *isRequired)
try {
  ensureConstructed();
  *isRequired = impl->requiresWritingCheckpoint() ? 1 : 0;
} catch (::precice::Error &e) {
  std::abort();
}

void precicef_requires_reading_checkpoint_(
    int *isRequired)
try {
  ensureConstructed();
  *isRequired = impl->requiresReadingCheckpoint() ? 1 : 0;
} catch (::precice::Error &e) {
  std::abort();
}

void precicef_requires_mesh_connectivity_for_(
    const char *meshName,
    int        *required,
    int         meshNameLength)
try {
  ensureConstructed();
  if (impl->requiresMeshConnectivityFor(strippedStringView(meshName, meshNameLength))) {
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
  ensureConstructed();
  impl->resetMesh(strippedStringView(meshName, meshNameLength));
} catch (::precice::Error &e) {
  std::abort();
}

void precicef_set_mesh_vertex_(
    const char   *meshName,
    const double *position,
    int          *vertexID,
    int           meshNameLength)
try {
  ensureConstructed();
  auto sv           = strippedStringView(meshName, meshNameLength);
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
  ensureConstructed();
  *meshSize = impl->getMeshVertexSize(strippedStringView(meshName, meshNameLength));
} catch (::precice::Error &e) {
  std::abort();
}

void precicef_set_mesh_vertices_(
    const char *meshName,
    const int  *size,
    double     *coordinates,
    int        *ids,
    int         meshNameLength)
try {
  ensureConstructed();
  auto sv           = strippedStringView(meshName, meshNameLength);
  auto positionSize = static_cast<unsigned long>(impl->getMeshDimensions(sv) * *size);
  impl->setMeshVertices(sv, {coordinates, positionSize}, {ids, static_cast<unsigned long>(*size)});
} catch (::precice::Error &e) {
  std::abort();
}

void precicef_set_mesh_edge_(
    const char *meshName,
    const int  *firstVertexID,
    const int  *secondVertexID,
    int         meshNameLength)
try {
  ensureConstructed();
  impl->setMeshEdge(strippedStringView(meshName, meshNameLength), *firstVertexID, *secondVertexID);
} catch (::precice::Error &e) {
  std::abort();
}

void precicef_set_mesh_edges_(
    const char *meshName,
    const int  *size,
    const int  *ids,
    int         meshNameLength)
try {
  ensureConstructed();
  auto idsSize = static_cast<unsigned long>(*size) * 2;
  impl->setMeshEdges(strippedStringView(meshName, meshNameLength), {ids, idsSize});
} catch (::precice::Error &e) {
  std::abort();
}

void precicef_set_mesh_triangle_(
    const char *meshName,
    const int  *firstVertexID,
    const int  *secondVertexID,
    const int  *thirdVertexID,
    int         meshNameLength)
try {
  ensureConstructed();
  impl->setMeshTriangle(strippedStringView(meshName, meshNameLength), *firstVertexID, *secondVertexID, *thirdVertexID);
} catch (::precice::Error &e) {
  std::abort();
}

void precicef_set_mesh_triangles_(
    const char *meshName,
    const int  *size,
    const int  *ids,
    int         meshNameLength)
try {
  ensureConstructed();
  auto idsSize = static_cast<unsigned long>(*size) * 3;
  impl->setMeshTriangles(strippedStringView(meshName, meshNameLength), {ids, idsSize});
} catch (::precice::Error &e) {
  std::abort();
}

void precicef_set_mesh_quad_(
    const char *meshName,
    const int  *firstVertexID,
    const int  *secondVertexID,
    const int  *thirdVertexID,
    const int  *fourthVertexID,
    int         meshNameLength)
try {
  ensureConstructed();
  impl->setMeshQuad(strippedStringView(meshName, meshNameLength), *firstVertexID, *secondVertexID, *thirdVertexID, *fourthVertexID);
} catch (::precice::Error &e) {
  std::abort();
}

void precicef_set_mesh_quads_(
    const char *meshName,
    const int  *size,
    const int  *ids,
    int         meshNameLength)
try {
  ensureConstructed();
  auto idsSize = static_cast<unsigned long>(*size) * 4;
  impl->setMeshQuads(strippedStringView(meshName, meshNameLength), {ids, idsSize});
} catch (::precice::Error &e) {
  std::abort();
}

void precicef_set_mesh_tetrahedron_(
    const char *meshName,
    const int  *firstVertexID,
    const int  *secondVertexID,
    const int  *thirdVertexID,
    const int  *fourthVertexID,
    int         meshNameLength)
try {
  ensureConstructed();
  impl->setMeshTetrahedron(strippedStringView(meshName, meshNameLength), *firstVertexID, *secondVertexID, *thirdVertexID, *fourthVertexID);
} catch (::precice::Error &e) {
  std::abort();
}

void precicef_set_mesh_tetrahedra_(
    const char *meshName,
    const int  *size,
    const int  *ids,
    int         meshNameLength)
try {
  ensureConstructed();
  auto idsSize = static_cast<unsigned long>(*size) * 4;
  impl->setMeshTetrahedra(strippedStringView(meshName, meshNameLength), {ids, idsSize});
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
  ensureConstructed();
  auto strippedMeshName = strippedStringView(meshName, meshNameLength);
  auto strippedDataName = strippedStringView(dataName, dataNameLength);
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
  ensureConstructed();
  auto strippedMeshName = strippedStringView(meshName, meshNameLength);
  auto strippedDataName = strippedStringView(dataName, dataNameLength);
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
  ensureConstructed();
  if (impl->requiresGradientDataFor(strippedStringView(meshName, meshNameLength), strippedStringView(dataName, dataNameLength))) {
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
  ensureConstructed();
  auto strippedMeshName   = strippedStringView(meshName, meshNameLength);
  auto strippedDataName   = strippedStringView(dataName, dataNameLength);
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
  ensureConstructed();
  auto strippedMeshName = strippedStringView(meshName, meshNameLength);
  auto strippedDataName = strippedStringView(dataName, dataNameLength);
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
  ensureConstructed();
  auto strippedMeshName = strippedStringView(meshName, meshNameLength);
  auto strippedDataName = strippedStringView(dataName, dataNameLength);
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
  ensureConstructed();
  auto sv     = strippedStringView(meshName, meshNameLength);
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
  ensureConstructed();
  auto sv              = strippedStringView(meshName, meshNameLength);
  auto coordinatesSize = static_cast<unsigned long>(impl->getMeshDimensions(sv) * *size);
  impl->getMeshVertexIDsAndCoordinates(sv, {ids, static_cast<unsigned long>(*size)}, {coordinates, coordinatesSize});
} catch (::precice::Error &e) {
  std::abort();
}

void precicef_start_profiling_section_(
    const char *sectionName,
    int         sectionNameLength)
try {
  ensureConstructed();
  auto sv = strippedStringView(sectionName, sectionNameLength);
  impl->startProfilingSection(sv);
} catch (::precice::Error &e) {
  std::abort();
}

void precicef_stop_last_profiling_section_()
try {
  ensureConstructed();
  impl->stopLastProfilingSection();
} catch (::precice::Error &e) {
  std::abort();
}

// Deprecated wrappers for backward compatibility

void precicef_set_vertex_(
    const char   *meshName,
    const double *coordinates,
    int          *id,
    int           meshNameLength)
{
  precicef_set_mesh_vertex_(meshName, coordinates, id, meshNameLength);
}

void precicef_set_edge_(
    const char *meshName,
    const int  *firstVertexID,
    const int  *secondVertexID,
    int         meshNameLength)
{
  precicef_set_mesh_edge_(meshName, firstVertexID, secondVertexID, meshNameLength);
}

void precicef_set_triangle_(
    const char *meshName,
    const int  *firstVertexID,
    const int  *secondVertexID,
    const int  *thirdVertexID,
    int         meshNameLength)
{
  precicef_set_mesh_triangle_(meshName, firstVertexID, secondVertexID, thirdVertexID, meshNameLength);
}

void precicef_set_quad_(
    const char *meshName,
    const int  *firstVertexID,
    const int  *secondVertexID,
    const int  *thirdVertexID,
    const int  *fourthVertexID,
    int         meshNameLength)
{
  precicef_set_mesh_quad_(meshName, firstVertexID, secondVertexID, thirdVertexID, fourthVertexID, meshNameLength);
}

void precicef_set_tetrahedron(
    const char *meshName,
    const int  *firstVertexID,
    const int  *secondVertexID,
    const int  *thirdVertexID,
    const int  *fourthVertexID,
    int         meshNameLength)
{
  precicef_set_mesh_tetrahedron_(meshName, firstVertexID, secondVertexID, thirdVertexID, fourthVertexID, meshNameLength);
}

void precicef_set_vertices_(
    const char *meshName,
    const int  *size,
    double     *coordinates,
    int        *ids,
    int         meshNameLength)
{
  precicef_set_mesh_vertices_(meshName, size, coordinates, ids, meshNameLength);
}

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif
#ifdef __clang__
#pragma clang diagnostic pop
#endif
