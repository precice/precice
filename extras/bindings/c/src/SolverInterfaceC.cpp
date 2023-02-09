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
  std::string stringAccessorName(participantName);
  std::string stringConfigFileName(configFileName);

  PRECICE_CHECK(impl == nullptr, errormsgCreate);
  impl.reset(new precice::SolverInterface(stringAccessorName,
                                          stringConfigFileName,
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
  std::string stringAccessorName(participantName);
  std::string stringConfigFileName(configFileName);

  PRECICE_CHECK(impl == nullptr, errormsgCreate);
  impl.reset(new precice::SolverInterface(stringAccessorName,
                                          stringConfigFileName,
                                          solverProcessIndex,
                                          solverProcessSize));
}

double precicec_initialize()
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  return impl->initialize();
}

double precicec_advance(double computedTimestepLength)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  return impl->advance(computedTimestepLength);
}

void precicec_finalize()
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->finalize();
  impl.reset();
}

int precicec_getDimensions()
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  return impl->getDimensions();
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
  std::string stringMeshName(meshName);
  if (impl->hasMesh(stringMeshName)) {
    return 1;
  }
  return 0;
}

int precicec_getMeshID(const char *meshName)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  std::string stringMeshName(meshName);
  return impl->getMeshID(stringMeshName);
}

int precicec_hasData(const char *dataName, int meshID)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  std::string stringDataName(dataName);
  return impl->hasData(stringDataName, meshID);
}

int precicec_getDataID(const char *dataName, int meshID)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  std::string stringDataName(dataName);
  return impl->getDataID(stringDataName, meshID);
}

int precicec_requiresMeshConnectivityFor(int meshID)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  if (impl->requiresMeshConnectivityFor(meshID)) {
    return 1;
  }
  return 0;
}

int precicec_setMeshVertex(
    int           meshID,
    const double *position)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  return impl->setMeshVertex(meshID, position);
}

void precicec_setMeshVertices(
    int           meshID,
    int           size,
    const double *positions,
    int *         ids)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->setMeshVertices(meshID, size, positions, ids);
}

int precicec_getMeshVertexSize(
    int meshID)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  return impl->getMeshVertexSize(meshID);
}

void precicec_setMeshEdge(
    int meshID,
    int firstVertexID,
    int secondVertexID)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->setMeshEdge(meshID, firstVertexID, secondVertexID);
}

void precicec_setMeshEdges(
    int        meshID,
    int        size,
    const int *vertices)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->setMeshEdges(meshID, size, vertices);
}

void precicec_setMeshTriangle(
    int meshID,
    int firstVertexID,
    int secondVertexID,
    int thirdVertexID)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->setMeshTriangle(meshID, firstVertexID, secondVertexID, thirdVertexID);
}

void precicec_setMeshTriangles(
    int        meshID,
    int        size,
    const int *vertices)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->setMeshTriangles(meshID, size, vertices);
}

void precicec_setMeshQuad(
    int meshID,
    int firstVertexID,
    int secondVertexID,
    int thirdVertexID,
    int fourthVertexID)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->setMeshQuad(meshID, firstVertexID, secondVertexID, thirdVertexID, fourthVertexID);
}

void precicec_setMeshQuads(
    int        meshID,
    int        size,
    const int *vertices)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->setMeshQuads(meshID, size, vertices);
}

void precicec_setMeshTetrahedron(
    int meshID,
    int firstVertexID,
    int secondVertexID,
    int thirdVertexID,
    int fourthVertexID)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->setMeshTetrahedron(meshID, firstVertexID, secondVertexID, thirdVertexID, fourthVertexID);
}

void precicec_setMeshTetrahedra(
    int        meshID,
    int        size,
    const int *vertices)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->setMeshTetrahedra(meshID, size, vertices);
}

void precicec_writeBlockVectorData(
    int           dataID,
    int           size,
    const int *   valueIndices,
    const double *values)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->writeBlockVectorData(dataID, size, valueIndices, values);
}

void precicec_writeVectorData(
    int           dataID,
    int           valueIndex,
    const double *dataValue)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->writeVectorData(dataID, valueIndex, dataValue);
}

void precicec_writeBlockScalarData(
    int           dataID,
    int           size,
    const int *   valueIndices,
    const double *values)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->writeBlockScalarData(dataID, size, valueIndices, values);
}

void precicec_writeScalarData(
    int    dataID,
    int    valueIndex,
    double dataValue)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->writeScalarData(dataID, valueIndex, dataValue);
}

void precicec_readBlockVectorData(
    int        dataID,
    int        size,
    const int *valueIndices,
    double *   values)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->readBlockVectorData(dataID, size, valueIndices, values);
}

void precicec_readVectorData(
    int     dataID,
    int     valueIndex,
    double *dataValue)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->readVectorData(dataID, valueIndex, dataValue);
}

void precicec_readBlockScalarData(
    int        dataID,
    int        size,
    const int *valueIndices,
    double *   values)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->readBlockScalarData(dataID, size, valueIndices, values);
}

void precicec_readScalarData(
    int     dataID,
    int     valueIndex,
    double *dataValue)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->readScalarData(dataID, valueIndex, *dataValue);
}

int precicec_requiresGradientDataFor(int dataID)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  if (impl->requiresGradientDataFor(dataID)) {
    return 1;
  }
  return 0;
}

void precicec_writeScalarGradientData(
    int           dataID,
    int           valueIndex,
    const double *gradientValues)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->writeScalarGradientData(dataID, valueIndex, gradientValues);
}

void precicec_writeBlockScalarGradientData(
    int           dataID,
    int           size,
    const int *   valueIndices,
    const double *gradientValues)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->writeBlockScalarGradientData(dataID, size, valueIndices, gradientValues);
}

void precicec_writeVectorGradientData(
    int           dataID,
    int           valueIndex,
    const double *gradientValues)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->writeVectorGradientData(dataID, valueIndex, gradientValues);
}

void precicec_writeBlockVectorGradientData(
    int           dataID,
    int           size,
    const int *   valueIndices,
    const double *gradientValues)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->writeBlockVectorGradientData(dataID, size, valueIndices, gradientValues);
}

const char *precicec_getVersionInformation()
{
  return precice::versionInformation;
}

void precicec_setMeshAccessRegion(
    const int     meshID,
    const double *boundingBox)
{
  impl->setMeshAccessRegion(meshID, boundingBox);
}

void precicec_getMeshVerticesAndIDs(
    const int meshID,
    const int size,
    int *     ids,
    double *  coordinates)
{
  impl->getMeshVerticesAndIDs(meshID, size, ids, coordinates);
}

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif
#ifdef __clang__
#pragma clang diagnostic pop
#endif
