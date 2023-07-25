extern "C" {
#include "precice/preciceC.h"
}
#include <memory>
#include <string>
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

static std::unique_ptr<precice::Participant> impl = nullptr;

static precice::logging::Logger _log("precicec");

static std::string errormsg       = "preCICE has not been created properly. Be sure to call \"precicec_createParticipant\" or \"precicec_createParticipant_withCommunicator\" before any other call to preCICE.";
static std::string errormsgCreate = "preCICE has been created already! Be sure to call either \"precicec_createParticipant\" or \"precicec_createParticipant_withCommunicator\" exactly once.";

void precicec_createParticipant_withCommunicator(
    const char *participantName,
    const char *configFileName,
    int         solverProcessIndex,
    int         solverProcessSize,
    void *      communicator)
{
  PRECICE_CHECK(impl == nullptr, errormsgCreate);
  impl.reset(new precice::Participant(participantName,
                                      configFileName,
                                      solverProcessIndex,
                                      solverProcessSize,
                                      communicator));
}

void precicec_createParticipant(
    const char *participantName,
    const char *configFileName,
    int         solverProcessIndex,
    int         solverProcessSize)
{
  PRECICE_CHECK(impl == nullptr, errormsgCreate);
  impl.reset(new precice::Participant(participantName,
                                      configFileName,
                                      solverProcessIndex,
                                      solverProcessSize));
}

void precicec_initialize()
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->initialize();
}

void precicec_advance(double computedTimeStepSize)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->advance(computedTimeStepSize);
}

void precicec_finalize()
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->finalize();
  impl.reset();
}

int precicec_getMeshDimensions(const char *meshName)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  return impl->getMeshDimensions(meshName);
}

int precicec_getDataDimensions(const char *meshName, const char *dataName)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  return impl->getDataDimensions(meshName, dataName);
}

int precicec_isCouplingOngoing()
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  if (impl->isCouplingOngoing()) {
    return 1;
  }
  return 0;
}

int precicec_isTimeWindowComplete()
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  if (impl->isTimeWindowComplete()) {
    return 1;
  }
  return 0;
}

double precicec_getMaxTimeStepSize()
{
  return impl->getMaxTimeStepSize();
}

int precicec_requiresInitialData()
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  return impl->requiresInitialData() ? 1 : 0;
}

int precicec_requiresWritingCheckpoint()
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  return impl->requiresWritingCheckpoint() ? 1 : 0;
}

int precicec_requiresReadingCheckpoint()
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  return impl->requiresReadingCheckpoint() ? 1 : 0;
}

int precicec_hasMesh(const char *meshName)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  if (impl->hasMesh(meshName)) {
    return 1;
  }
  return 0;
}

int precicec_hasData(const char *meshName, const char *dataName)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  return impl->hasData(meshName, dataName);
}

int precicec_requiresMeshConnectivityFor(const char *meshName)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  if (impl->requiresMeshConnectivityFor(meshName)) {
    return 1;
  }
  return 0;
}

int precicec_setMeshVertex(
    const char *  meshName,
    const double *position)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  auto size = static_cast<long unsigned>(impl->getMeshDimensions(meshName));
  return impl->setMeshVertex(meshName, {position, size});
}

void precicec_setMeshVertices(
    const char *  meshName,
    int           size,
    const double *positions,
    int *         ids)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  auto idsSize = static_cast<long unsigned>(size);
  auto posSize = static_cast<long unsigned>(impl->getMeshDimensions(meshName) * size);
  impl->setMeshVertices(meshName, {positions, posSize}, {ids, idsSize});
}

int precicec_getMeshVertexSize(
    const char *meshName)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  return impl->getMeshVertexSize(meshName);
}

void precicec_setMeshEdge(
    const char *meshName,
    int         firstVertexID,
    int         secondVertexID)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->setMeshEdge(meshName, firstVertexID, secondVertexID);
}

void precicec_setMeshEdges(
    const char *meshName,
    int         size,
    const int * vertices)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  auto verticesSize = static_cast<long unsigned>(size) * 2;
  impl->setMeshEdges(meshName, {vertices, verticesSize});
}

void precicec_setMeshTriangle(
    const char *meshName,
    int         firstVertexID,
    int         secondVertexID,
    int         thirdVertexID)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->setMeshTriangle(meshName, firstVertexID, secondVertexID, thirdVertexID);
}

void precicec_setMeshTriangles(
    const char *meshName,
    int         size,
    const int * vertices)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  auto verticesSize = static_cast<long unsigned>(size) * 3;
  impl->setMeshTriangles(meshName, {vertices, verticesSize});
}

void precicec_setMeshQuad(
    const char *meshName,
    int         firstVertexID,
    int         secondVertexID,
    int         thirdVertexID,
    int         fourthVertexID)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->setMeshQuad(meshName, firstVertexID, secondVertexID, thirdVertexID, fourthVertexID);
}

void precicec_setMeshQuads(
    const char *meshName,
    int         size,
    const int * vertices)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  auto verticesSize = static_cast<long unsigned>(size) * 4;
  impl->setMeshQuads(meshName, {vertices, verticesSize});
}

void precicec_setMeshTetrahedron(
    const char *meshName,
    int         firstVertexID,
    int         secondVertexID,
    int         thirdVertexID,
    int         fourthVertexID)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->setMeshTetrahedron(meshName, firstVertexID, secondVertexID, thirdVertexID, fourthVertexID);
}

void precicec_setMeshTetrahedra(
    const char *meshName,
    int         size,
    const int * vertices)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  auto verticesSize = static_cast<long unsigned>(size) * 4;
  impl->setMeshTetrahedra(meshName, {vertices, verticesSize});
}

void precicec_writeData(
    const char *  meshName,
    const char *  dataName,
    int           size,
    const int *   valueIndices,
    const double *values)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  auto dataSize = size * impl->getDataDimensions(meshName, dataName);
  impl->writeData(meshName, dataName, {valueIndices, static_cast<unsigned long>(size)}, {values, static_cast<unsigned long>(dataSize)});
}

void precicec_readData(
    const char *meshName,
    const char *dataName,
    int         size,
    const int * valueIndices,
    double      relativeReadTime,
    double *    values)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  auto dataSize = size * impl->getDataDimensions(meshName, dataName);
  impl->readData(meshName, dataName, {valueIndices, static_cast<unsigned long>(size)}, relativeReadTime, {values, static_cast<unsigned long>(dataSize)});
}

int precicec_requiresGradientDataFor(const char *meshName,
                                     const char *dataName)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  if (impl->requiresGradientDataFor(meshName, dataName)) {
    return 1;
  }
  return 0;
}

void precicec_writeGradientData(
    const char *  meshName,
    const char *  dataName,
    int           size,
    const int *   valueIndices,
    const double *gradients)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  auto gradientComponents = impl->getDataDimensions(meshName, dataName) * impl->getMeshDimensions(meshName);
  auto gradientSize       = size * gradientComponents;
  impl->writeGradientData(meshName, dataName, {valueIndices, static_cast<unsigned long>(size)}, {gradients, static_cast<unsigned long>(gradientSize)});
}

const char *precicec_getVersionInformation()
{
  return precice::versionInformation;
}

void precicec_setMeshAccessRegion(
    const char *  meshName,
    const double *boundingBox)
{
  auto bbSize = static_cast<long unsigned>(impl->getMeshDimensions(meshName)) * 2;
  impl->setMeshAccessRegion(meshName, {boundingBox, bbSize});
}

void precicec_getMeshVertexIDsAndCoordinates(
    const char *meshName,
    const int   size,
    int *       ids,
    double *    coordinates)
{
  auto coordinatesSize = static_cast<long unsigned>(impl->getMeshDimensions(meshName) * size);
  impl->getMeshVertexIDsAndCoordinates(meshName, {ids, static_cast<unsigned long>(size)}, {coordinates, coordinatesSize});
}

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif
#ifdef __clang__
#pragma clang diagnostic pop
#endif
