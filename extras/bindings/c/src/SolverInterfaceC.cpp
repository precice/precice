extern "C" {
#include "precice/SolverInterfaceC.h"
}
#include <memory>
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

static std::unique_ptr<precice::SolverInterface> impl = nullptr;

static precice::logging::Logger _log("SolverInterfaceC");

static std::string errormsg       = "preCICE has not been created properly. Be sure to call \"precicec_createSolverInterface\" or \"precicec_createSolverInterface_withCommunicator\" before any other call to preCICE.";
static std::string errormsgCreate = "preCICE has been created already! Be sure to call either \"precicec_createSolverInterface\" or \"precicec_createSolverInterface_withCommunicator\" exactly once.";

void precicec_createSolverInterface_withCommunicator(
    const char *participantName,
    const char *configFileName,
    int         solverProcessIndex,
    int         solverProcessSize,
    void *      communicator)
{
  PRECICE_CHECK(impl == nullptr, errormsgCreate);
  impl.reset(new precice::SolverInterface(participantName,
                                          configFileName,
                                          solverProcessIndex,
                                          solverProcessSize,
                                          communicator));
}

void precicec_createSolverInterface(
    const char *participantName,
    const char *configFileName,
    int         solverProcessIndex,
    int         solverProcessSize)
{
  PRECICE_CHECK(impl == nullptr, errormsgCreate);
  impl.reset(new precice::SolverInterface(participantName,
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
  return impl->setMeshVertex(meshName, position);
}

void precicec_setMeshVertices(
    const char *  meshName,
    int           size,
    const double *positions,
    int *         ids)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->setMeshVertices(meshName, size, positions, ids);
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
  impl->setMeshEdges(meshName, size, vertices);
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
  impl->setMeshTriangles(meshName, size, vertices);
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
  impl->setMeshQuads(meshName, size, vertices);
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
  impl->setMeshTetrahedra(meshName, size, vertices);
}

void precicec_writeBlockVectorData(
    const char *  meshName,
    const char *  dataName,
    int           size,
    const int *   valueIndices,
    const double *values)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->writeBlockVectorData(meshName, dataName, size, valueIndices, values);
}

void precicec_writeVectorData(
    const char *  meshName,
    const char *  dataName,
    int           valueIndex,
    const double *dataValue)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->writeVectorData(meshName, dataName, valueIndex, dataValue);
}

void precicec_writeBlockScalarData(
    const char *  meshName,
    const char *  dataName,
    int           size,
    const int *   valueIndices,
    const double *values)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->writeBlockScalarData(meshName, dataName, size, valueIndices, values);
}

void precicec_writeScalarData(
    const char *meshName,
    const char *dataName,
    int         valueIndex,
    double      dataValue)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->writeScalarData(meshName, dataName, valueIndex, dataValue);
}

void precicec_readBlockVectorData(
    const char *meshName,
    const char *dataName,
    int         size,
    const int * valueIndices,
    double      relativeReadTime,
    double *    values)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->readBlockVectorData(meshName, dataName, size, valueIndices, relativeReadTime, values);
}

void precicec_readVectorData(
    const char *meshName,
    const char *dataName,
    int         valueIndex,
    double      relativeReadTime,
    double *    dataValue)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->readVectorData(meshName, dataName, valueIndex, relativeReadTime, dataValue);
}

void precicec_readBlockScalarData(
    const char *meshName,
    const char *dataName,
    int         size,
    const int * valueIndices,
    double      relativeReadTime,
    double *    values)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->readBlockScalarData(meshName, dataName, size, valueIndices, relativeReadTime, values);
}

void precicec_readScalarData(
    const char *meshName,
    const char *dataName,
    int         valueIndex,
    double      relativeReadTime,
    double *    dataValue)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->readScalarData(meshName, dataName, valueIndex, relativeReadTime, *dataValue);
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

void precicec_writeScalarGradientData(
    const char *  meshName,
    const char *  dataName,
    int           valueIndex,
    const double *gradientValues)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->writeScalarGradientData(meshName, dataName, valueIndex, gradientValues);
}

void precicec_writeBlockScalarGradientData(
    const char *  meshName,
    const char *  dataName,
    int           size,
    const int *   valueIndices,
    const double *gradientValues)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->writeBlockScalarGradientData(meshName, dataName, size, valueIndices, gradientValues);
}

void precicec_writeVectorGradientData(
    const char *  meshName,
    const char *  dataName,
    int           valueIndex,
    const double *gradientValues)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->writeVectorGradientData(meshName, dataName, valueIndex, gradientValues);
}

void precicec_writeBlockVectorGradientData(
    const char *  meshName,
    const char *  dataName,
    int           size,
    const int *   valueIndices,
    const double *gradientValues)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->writeBlockVectorGradientData(meshName, dataName, size, valueIndices, gradientValues);
}

const char *precicec_getVersionInformation()
{
  return precice::versionInformation;
}

void precicec_setMeshAccessRegion(
    const char *  meshName,
    const double *boundingBox)
{
  impl->setMeshAccessRegion(meshName, boundingBox);
}

void precicec_getMeshVerticesAndIDs(
    const char *meshName,
    const int   size,
    int *       ids,
    double *    coordinates)
{
  impl->getMeshVerticesAndIDs(meshName, size, ids, coordinates);
}

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif
#ifdef __clang__
#pragma clang diagnostic pop
#endif
