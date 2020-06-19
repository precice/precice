extern "C" {
#include "precice/SolverInterfaceC.h"
}
#include <memory>
#include <string>
#include "precice/SolverInterface.hpp"
#include "precice/impl/versions.hpp"
#include "utils/assertion.hpp"

static std::unique_ptr<precice::SolverInterface> interface;

void precicec_createSolverInterface(
    const char *participantName,
    const char *configFileName,
    int         solverProcessIndex,
    int         solverProcessSize)
{
  std::string stringAccessorName(participantName);
  std::string stringConfigFileName(configFileName);

  PRECICE_ASSERT(interface == nullptr);
  interface = std::unique_ptr<precice::SolverInterface>( new precice::SolverInterface(stringAccessorName,
                                                                                     stringConfigFileName,
                                                                                     solverProcessIndex,
                                                                                     solverProcessSize )
                                                        );
}

void precicec_createSolverInterface_withCommunicator(
    const char *participantName,
    const char *configFileName,
    int         solverProcessIndex,
    int         solverProcessSize,
    void *      communicator)
{
  std::string stringAccessorName(participantName);
  std::string stringConfigFileName(configFileName);

  PRECICE_ASSERT(interface == nullptr);
  interface = std::unique_ptr<precice::SolverInterface>( new precice::SolverInterface(stringAccessorName,
                                                                                     stringConfigFileName,
                                                                                     solverProcessIndex,
                                                                                     solverProcessSize,
                                                                                     communicator)
                                                        );
}

double precicec_initialize()
{
  PRECICE_ASSERT(interface != nullptr);
  return interface->initialize();
}

void precicec_initialize_data()
{
  PRECICE_ASSERT(interface != nullptr);
  interface->initializeData();
}

double precicec_advance(double computedTimestepLength)
{
  PRECICE_ASSERT(interface != nullptr);
  return interface->advance(computedTimestepLength);
}

void precicec_finalize()
{
  if (interface != nullptr)
  {
    interface->finalize();
  }
}

int precicec_getDimensions()
{
  PRECICE_ASSERT(interface != nullptr);
  return interface->getDimensions();
}

int precicec_isCouplingOngoing()
{
  PRECICE_ASSERT(interface != nullptr);
  if (interface->isCouplingOngoing()) {
    return 1;
  }
  return 0;
}

int precicec_isTimeWindowComplete()
{
  PRECICE_ASSERT(interface != nullptr);
  if (interface->isTimeWindowComplete()) {
    return 1;
  }
  return 0;
}

int precicec_hasToEvaluateSurrogateModel()
{
  PRECICE_ASSERT(interface != nullptr);
  if (interface->hasToEvaluateSurrogateModel()) {
    return 1;
  }
  return 0;
}

int precicec_hasToEvaluateFineModel()
{
  PRECICE_ASSERT(interface != nullptr);
  if (interface->hasToEvaluateFineModel()) {
    return 1;
  }
  return 0;
}

int precicec_isReadDataAvailable()
{
  PRECICE_ASSERT(interface != nullptr);
  if (interface->isReadDataAvailable()) {
    return 1;
  }
  return 0;
}

int precicec_isWriteDataRequired(double computedTimestepLength)
{
  PRECICE_ASSERT(interface != nullptr);
  if (interface->isWriteDataRequired(computedTimestepLength)) {
    return 1;
  }
  return 0;
}

int precicec_isActionRequired(const char *action)
{
  PRECICE_ASSERT(interface != nullptr);
  PRECICE_ASSERT(action != nullptr);
  if (interface->isActionRequired(std::string(action))) {
    return 1;
  }
  return 0;
}

void precicec_markActionFulfilled(const char *action)
{
  PRECICE_ASSERT(interface != nullptr);
  PRECICE_ASSERT(action != nullptr);
  interface->markActionFulfilled(std::string(action));
}

int precicec_hasMesh(const char *meshName)
{
  PRECICE_ASSERT(interface != nullptr);
  std::string stringMeshName(meshName);
  if (interface->hasMesh(stringMeshName)) {
    return 1;
  }
  return 0;
}

int precicec_getMeshID(const char *meshName)
{
  PRECICE_ASSERT(interface != nullptr);
  std::string stringMeshName(meshName);
  return interface->getMeshID(stringMeshName);
}

int precicec_hasData(const char *dataName, int meshID)
{
  PRECICE_ASSERT(interface != nullptr);
  std::string stringDataName(dataName);
  return interface->hasData(stringDataName, meshID);
}

int precicec_getDataID(const char *dataName, int meshID)
{
  PRECICE_ASSERT(interface != nullptr);
  std::string stringDataName(dataName);
  return interface->getDataID(stringDataName, meshID);
}

int precicec_setMeshVertex(
    int           meshID,
    const double *position)
{
  PRECICE_ASSERT(interface != nullptr);
  return interface->setMeshVertex(meshID, position);
}

void precicec_getMeshVertices(
    int        meshID,
    int        size,
    const int *ids,
    double *   positions)
{
  PRECICE_ASSERT(interface != nullptr);
  interface->getMeshVertices(meshID, size, ids, positions);
}

void precicec_setMeshVertices(
    int           meshID,
    int           size,
    const double *positions,
    int *         ids)
{
  PRECICE_ASSERT(interface != nullptr);
  interface->setMeshVertices(meshID, size, positions, ids);
}

int precicec_getMeshVertexSize(
    int meshID)
{
  PRECICE_ASSERT(interface != nullptr);
  return interface->getMeshVertexSize(meshID);
}

void precicec_getMeshVertexIDsFromPositions(
    int           meshID,
    int           size,
    const double *positions,
    int *         ids)
{
  PRECICE_ASSERT(interface != nullptr);
  interface->getMeshVertexIDsFromPositions(meshID, size, positions, ids);
}

int precicec_setMeshEdge(
    int meshID,
    int firstVertexID,
    int secondVertexID)
{
  PRECICE_ASSERT(interface != nullptr);
  return interface->setMeshEdge(meshID, firstVertexID, secondVertexID);
}

void precicec_setMeshTriangle(
    int meshID,
    int firstEdgeID,
    int secondEdgeID,
    int thirdEdgeID)
{
  PRECICE_ASSERT(interface != nullptr);
  interface->setMeshTriangle(meshID, firstEdgeID, secondEdgeID, thirdEdgeID);
}

void precicec_setMeshTriangleWithEdges(
    int meshID,
    int firstVertexID,
    int secondVertexID,
    int thirdVertexID)
{
  PRECICE_ASSERT(interface != nullptr);
  interface->setMeshTriangleWithEdges(meshID, firstVertexID, secondVertexID, thirdVertexID);
}

void precicec_setMeshQuad(
    int meshID,
    int firstEdgeID,
    int secondEdgeID,
    int thirdEdgeID,
    int fourthEdgeID)
{
  PRECICE_ASSERT(interface != nullptr);
  interface->setMeshQuad(meshID, firstEdgeID, secondEdgeID, thirdEdgeID, fourthEdgeID);
}

void precicec_setMeshQuadWithEdges(
    int meshID,
    int firstVertexID,
    int secondVertexID,
    int thirdVertexID,
    int fourthVertexID)
{
  PRECICE_ASSERT(interface != nullptr);
  interface->setMeshQuadWithEdges(meshID, firstVertexID, secondVertexID, thirdVertexID, fourthVertexID);
}

void precicec_writeBlockVectorData(
    int           dataID,
    int           size,
    const int *   valueIndices,
    const double *values)
{
  PRECICE_ASSERT(interface != nullptr);
  interface->writeBlockVectorData(dataID, size, valueIndices, values);
}

void precicec_writeVectorData(
    int           dataID,
    int           valueIndex,
    const double *dataValue)
{
  PRECICE_ASSERT(interface != nullptr);
  interface->writeVectorData(dataID, valueIndex, dataValue);
}

void precicec_writeBlockScalarData(
    int           dataID,
    int           size,
    const int *   valueIndices,
    const double *values)
{
  PRECICE_ASSERT(interface != nullptr);
  interface->writeBlockScalarData(dataID, size, valueIndices, values);
}

void precicec_writeScalarData(
    int    dataID,
    int    valueIndex,
    double dataValue)
{
  PRECICE_ASSERT(interface != nullptr);
  interface->writeScalarData(dataID, valueIndex, dataValue);
}

void precicec_readBlockVectorData(
    int        dataID,
    int        size,
    const int *valueIndices,
    double *   values)
{
  PRECICE_ASSERT(interface != nullptr);
  interface->readBlockVectorData(dataID, size, valueIndices, values);
}

void precicec_readVectorData(
    int     dataID,
    int     valueIndex,
    double *dataValue)
{
  PRECICE_ASSERT(interface != nullptr);
  interface->readVectorData(dataID, valueIndex, dataValue);
}

void precicec_readBlockScalarData(
    int        dataID,
    int        size,
    const int *valueIndices,
    double *   values)
{
  PRECICE_ASSERT(interface != nullptr);
  interface->readBlockScalarData(dataID, size, valueIndices, values);
}

void precicec_readScalarData(
    int     dataID,
    int     valueIndex,
    double *dataValue)
{
  PRECICE_ASSERT(interface != nullptr);
  interface->readScalarData(dataID, valueIndex, *dataValue);
}

const char *precicec_getVersionInformation()
{
  return precice::versionInformation;
}

void precicec_mapWriteDataFrom(int fromMeshID)
{
  PRECICE_ASSERT(interface != nullptr);
  interface->mapWriteDataFrom(fromMeshID);
}

void precicec_mapReadDataTo(int toMeshID)
{
  PRECICE_ASSERT(interface != nullptr);
  interface->mapReadDataTo(toMeshID);
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
