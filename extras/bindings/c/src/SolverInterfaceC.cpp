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

int precicec_hasData(const char *mesh, const char *data)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  return impl->hasData(mesh, data);
}

int precicec_requiresMeshConnectivityFor(const char *mesh)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  if (impl->requiresMeshConnectivityFor(mesh)) {
    return 1;
  }
  return 0;
}

int precicec_setMeshVertex(
    const char *  mesh,
    const double *position)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  return impl->setMeshVertex(mesh, position);
}

void precicec_setMeshVertices(
    const char *  mesh,
    int           size,
    const double *positions,
    int *         ids)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->setMeshVertices(mesh, size, positions, ids);
}

int precicec_getMeshVertexSize(
    const char *mesh)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  return impl->getMeshVertexSize(mesh);
}

void precicec_setMeshEdge(
    const char *mesh,
    int         firstVertexID,
    int         secondVertexID)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->setMeshEdge(mesh, firstVertexID, secondVertexID);
}

void precicec_setMeshEdges(
    const char *mesh,
    int         size,
    const int * vertices)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->setMeshEdges(mesh, size, vertices);
}

void precicec_setMeshTriangle(
    const char *mesh,
    int         firstVertexID,
    int         secondVertexID,
    int         thirdVertexID)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->setMeshTriangle(mesh, firstVertexID, secondVertexID, thirdVertexID);
}

void precicec_setMeshTriangles(
    const char *mesh,
    int         size,
    const int * vertices)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->setMeshTriangles(mesh, size, vertices);
}

void precicec_setMeshQuad(
    const char *mesh,
    int         firstVertexID,
    int         secondVertexID,
    int         thirdVertexID,
    int         fourthVertexID)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->setMeshQuad(mesh, firstVertexID, secondVertexID, thirdVertexID, fourthVertexID);
}

void precicec_setMeshQuads(
    const char *mesh,
    int         size,
    const int * vertices)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->setMeshQuads(mesh, size, vertices);
}

void precicec_setMeshTetrahedron(
    const char *mesh,
    int         firstVertexID,
    int         secondVertexID,
    int         thirdVertexID,
    int         fourthVertexID)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->setMeshTetrahedron(mesh, firstVertexID, secondVertexID, thirdVertexID, fourthVertexID);
}

void precicec_setMeshTetrahedra(
    const char *mesh,
    int         size,
    const int * vertices)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->setMeshTetrahedra(mesh, size, vertices);
}

void precicec_writeBlockVectorData(
    const char *  mesh,
    const char *  data,
    int           size,
    const int *   valueIndices,
    const double *values)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->writeBlockVectorData(mesh, data, size, valueIndices, values);
}

void precicec_writeVectorData(
    const char *  mesh,
    const char *  data,
    int           valueIndex,
    const double *dataValue)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->writeVectorData(mesh, data, valueIndex, dataValue);
}

void precicec_writeBlockScalarData(
    const char *  mesh,
    const char *  data,
    int           size,
    const int *   valueIndices,
    const double *values)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->writeBlockScalarData(mesh, data, size, valueIndices, values);
}

void precicec_writeScalarData(
    const char *mesh,
    const char *data,
    int         valueIndex,
    double      dataValue)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->writeScalarData(mesh, data, valueIndex, dataValue);
}

void precicec_readBlockVectorData(
    const char *mesh,
    const char *data,
    int         size,
    const int * valueIndices,
    double *    values)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->readBlockVectorData(mesh, data, size, valueIndices, values);
}

void precicec_readVectorData(
    const char *mesh,
    const char *data,
    int         valueIndex,
    double *    dataValue)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->readVectorData(mesh, data, valueIndex, dataValue);
}

void precicec_readBlockScalarData(
    const char *mesh,
    const char *data,
    int         size,
    const int * valueIndices,
    double *    values)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->readBlockScalarData(mesh, data, size, valueIndices, values);
}

void precicec_readScalarData(
    const char *mesh,
    const char *data,
    int         valueIndex,
    double *    dataValue)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->readScalarData(mesh, data, valueIndex, *dataValue);
}

int precicec_requiresGradientDataFor(const char *mesh,
                                     const char *data)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  if (impl->requiresGradientDataFor(mesh, data)) {
    return 1;
  }
  return 0;
}

void precicec_writeScalarGradientData(
    const char *  mesh,
    const char *  data,
    int           valueIndex,
    const double *gradientValues)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->writeScalarGradientData(mesh, data, valueIndex, gradientValues);
}

void precicec_writeBlockScalarGradientData(
    const char *  mesh,
    const char *  data,
    int           size,
    const int *   valueIndices,
    const double *gradientValues)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->writeBlockScalarGradientData(mesh, data, size, valueIndices, gradientValues);
}

void precicec_writeVectorGradientData(
    const char *  mesh,
    const char *  data,
    int           valueIndex,
    const double *gradientValues)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->writeVectorGradientData(mesh, data, valueIndex, gradientValues);
}

void precicec_writeBlockVectorGradientData(
    const char *  mesh,
    const char *  data,
    int           size,
    const int *   valueIndices,
    const double *gradientValues)
{
  PRECICE_CHECK(impl != nullptr, errormsg);
  impl->writeBlockVectorGradientData(mesh, data, size, valueIndices, gradientValues);
}

const char *precicec_getVersionInformation()
{
  return precice::versionInformation;
}

void precicec_setMeshAccessRegion(
    const char *  mesh,
    const double *boundingBox)
{
  impl->setMeshAccessRegion(mesh, boundingBox);
}

void precicec_getMeshVerticesAndIDs(
    const char *mesh,
    const int   size,
    int *       ids,
    double *    coordinates)
{
  impl->getMeshVerticesAndIDs(mesh, size, ids, coordinates);
}

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif
#ifdef __clang__
#pragma clang diagnostic pop
#endif
