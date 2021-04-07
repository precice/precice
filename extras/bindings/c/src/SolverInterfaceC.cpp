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
  precicec_createSolverInterface_withCommunicator(participantName,
                                                  configFileName,
                                                  solverProcessIndex,
                                                  solverProcessSize,
                                                  nullptr);
}

double precicec_initialize()
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  return impl->initialize();
}

void precicec_initialize_data()
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->initializeData();
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

int precicec_hasToEvaluateSurrogateModel()
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  if (impl->hasToEvaluateSurrogateModel()) {
    return 1;
  }
  return 0;
}

int precicec_hasToEvaluateFineModel()
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  if (impl->hasToEvaluateFineModel()) {
    return 1;
  }
  return 0;
}

int precicec_isReadDataAvailable()
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  if (impl->isReadDataAvailable()) {
    return 1;
  }
  return 0;
}

int precicec_isWriteDataRequired(double computedTimestepLength)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  if (impl->isWriteDataRequired(computedTimestepLength)) {
    return 1;
  }
  return 0;
}

int precicec_isActionRequired(const char *action)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  PRECICE_ASSERT(action != nullptr);
  if (impl->isActionRequired(std::string(action))) {
    return 1;
  }
  return 0;
}

void precicec_markActionFulfilled(const char *action)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  PRECICE_ASSERT(action != nullptr);
  impl->markActionFulfilled(std::string(action));
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

int precicec_setMeshVertex(
    int           meshID,
    const double *position)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  return impl->setMeshVertex(meshID, position);
}

void precicec_getMeshVertices(
    int        meshID,
    int        size,
    const int *ids,
    double *   positions)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->getMeshVertices(meshID, size, ids, positions);
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

void precicec_getMeshVertexIDsFromPositions(
    int           meshID,
    int           size,
    const double *positions,
    int *         ids)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->getMeshVertexIDsFromPositions(meshID, size, positions, ids);
}

int precicec_setMeshEdge(
    int meshID,
    int firstVertexID,
    int secondVertexID)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  return impl->setMeshEdge(meshID, firstVertexID, secondVertexID);
}

void precicec_setMeshTriangle(
    int meshID,
    int firstEdgeID,
    int secondEdgeID,
    int thirdEdgeID)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->setMeshTriangle(meshID, firstEdgeID, secondEdgeID, thirdEdgeID);
}

void precicec_setMeshTriangleWithEdges(
    int meshID,
    int firstVertexID,
    int secondVertexID,
    int thirdVertexID)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->setMeshTriangleWithEdges(meshID, firstVertexID, secondVertexID, thirdVertexID);
}

void precicec_setMeshQuad(
    int meshID,
    int firstEdgeID,
    int secondEdgeID,
    int thirdEdgeID,
    int fourthEdgeID)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->setMeshQuad(meshID, firstEdgeID, secondEdgeID, thirdEdgeID, fourthEdgeID);
}

void precicec_setMeshQuadWithEdges(
    int meshID,
    int firstVertexID,
    int secondVertexID,
    int thirdVertexID,
    int fourthVertexID)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->setMeshQuadWithEdges(meshID, firstVertexID, secondVertexID, thirdVertexID, fourthVertexID);
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

const char *precicec_getVersionInformation()
{
  return precice::versionInformation;
}

void precicec_mapWriteDataFrom(int fromMeshID)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->mapWriteDataFrom(fromMeshID);
}

void precicec_mapReadDataTo(int toMeshID)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->mapReadDataTo(toMeshID);
}

const char *precicec_actionWriteInitialData()
{
  return precice::constants::actionWriteInitialData().c_str();
}

const char *precicec_actionWriteIterationCheckpoint()
{
  return precice::constants::actionWriteIterationCheckpoint().c_str();
}

const char *precicec_actionReadIterationCheckpoint()
{
  return precice::constants::actionReadIterationCheckpoint().c_str();
}

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif
#ifdef __clang__
#pragma clang diagnostic pop
#endif
