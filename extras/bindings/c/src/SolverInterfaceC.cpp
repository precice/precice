extern "C" {
#include "precice/SolverInterfaceC.h"
}
#include <memory>
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

static std::unique_ptr<precice::SolverInterface> impl = nullptr;

static precice::logging::Logger _log("SolverInterfaceC");

static std::string errormsg       = "preCICE has not been created properly. Be sure to call \"precicec_createSolverInterface\" or \"precicec_createSolverInterface_withCommunicator\" before any other call to preCICE.";
static std::string errormsgCreate = "preCICE has been created already! Be sure to call either \"precicec_createSolverInterface\" or \"precicec_createSolverInterface_withCommunicator\" exactly once.";
using Err                         = ::precice::APIError;

void precicec_createSolverInterface_withCommunicator(
    const char *participantName,
    const char *configFileName,
    int         solverProcessIndex,
    int         solverProcessSize,
    void *      communicator)
try {
  std::string stringAccessorName(participantName);
  std::string stringConfigFileName(configFileName);

  PRECICE_CHECK(impl == nullptr, Err, errormsgCreate);
  impl.reset(new precice::SolverInterface(stringAccessorName,
                                          stringConfigFileName,
                                          solverProcessIndex,
                                          solverProcessSize,
                                          communicator));
} catch (...) {
  std::abort();
}

void precicec_createSolverInterface(
    const char *participantName,
    const char *configFileName,
    int         solverProcessIndex,
    int         solverProcessSize)
try {
  std::string stringAccessorName(participantName);
  std::string stringConfigFileName(configFileName);

  PRECICE_CHECK(impl == nullptr, Err, errormsgCreate);
  impl.reset(new precice::SolverInterface(stringAccessorName,
                                          stringConfigFileName,
                                          solverProcessIndex,
                                          solverProcessSize));
} catch (...) {
  std::abort();
}

double precicec_initialize()
try {
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  return impl->initialize();
} catch (...) {
  std::abort();
}

double precicec_advance(double computedTimestepLength)
try {
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  return impl->advance(computedTimestepLength);
} catch (...) {
  std::abort();
}

void precicec_finalize()
try {
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  impl->finalize();
  impl.reset();
} catch (...) {
  std::abort();
}

int precicec_getDimensions()
try {
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  return impl->getDimensions();
} catch (...) {
  std::abort();
}

int precicec_isCouplingOngoing()
try {
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  if (impl->isCouplingOngoing()) {
    return 1;
  }
  return 0;
} catch (...) {
  std::abort();
}

int precicec_isTimeWindowComplete()
try {
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  if (impl->isTimeWindowComplete()) {
    return 1;
  }
  return 0;
} catch (...) {
  std::abort();
}

int precicec_requiresInitialData()
try {
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  return impl->requiresInitialData() ? 1 : 0;
} catch (...) {
  std::abort();
}

int precicec_requiresWritingCheckpoint()
try {
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  return impl->requiresWritingCheckpoint() ? 1 : 0;
} catch (...) {
  std::abort();
}

int precicec_requiresReadingCheckpoint()
try {
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  return impl->requiresReadingCheckpoint() ? 1 : 0;
} catch (...) {
  std::abort();
}

int precicec_hasMesh(const char *meshName)
try {
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  std::string stringMeshName(meshName);
  if (impl->hasMesh(stringMeshName)) {
    return 1;
  }
  return 0;
} catch (...) {
  std::abort();
}

int precicec_getMeshID(const char *meshName)
try {
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  std::string stringMeshName(meshName);
  return impl->getMeshID(stringMeshName);
} catch (...) {
  std::abort();
}

int precicec_hasData(const char *dataName, int meshID)
try {
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  std::string stringDataName(dataName);
  return impl->hasData(stringDataName, meshID);
} catch (...) {
  std::abort();
}

int precicec_getDataID(const char *dataName, int meshID)
try {
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  std::string stringDataName(dataName);
  return impl->getDataID(stringDataName, meshID);
} catch (...) {
  std::abort();
}

int precicec_requiresMeshConnectivityFor(int meshID)
try {
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  if (impl->requiresMeshConnectivityFor(meshID)) {
    return 1;
  }
  return 0;
} catch (...) {
  std::abort();
}

int precicec_setMeshVertex(
    int           meshID,
    const double *position)
try {
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  return impl->setMeshVertex(meshID, position);
} catch (...) {
  std::abort();
}

void precicec_setMeshVertices(
    int           meshID,
    int           size,
    const double *positions,
    int *         ids)
try {
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  impl->setMeshVertices(meshID, size, positions, ids);
} catch (...) {
  std::abort();
}

int precicec_getMeshVertexSize(
    int meshID)
try {
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  return impl->getMeshVertexSize(meshID);
} catch (...) {
  std::abort();
}

void precicec_setMeshEdge(
    int meshID,
    int firstVertexID,
    int secondVertexID)
try {
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  impl->setMeshEdge(meshID, firstVertexID, secondVertexID);
} catch (...) {
  std::abort();
}

void precicec_setMeshEdges(
    int        meshID,
    int        size,
    const int *vertices)
try {
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  impl->setMeshEdges(meshID, size, vertices);
} catch (...) {
  std::abort();
}

void precicec_setMeshTriangle(
    int meshID,
    int firstVertexID,
    int secondVertexID,
    int thirdVertexID)
try {
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  impl->setMeshTriangle(meshID, firstVertexID, secondVertexID, thirdVertexID);
} catch (...) {
  std::abort();
}

void precicec_setMeshTriangles(
    int        meshID,
    int        size,
    const int *vertices)
try {
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  impl->setMeshTriangles(meshID, size, vertices);
} catch (...) {
  std::abort();
}

void precicec_setMeshQuad(
    int meshID,
    int firstVertexID,
    int secondVertexID,
    int thirdVertexID,
    int fourthVertexID)
try {
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  impl->setMeshQuad(meshID, firstVertexID, secondVertexID, thirdVertexID, fourthVertexID);
} catch (...) {
  std::abort();
}

void precicec_setMeshQuads(
    int        meshID,
    int        size,
    const int *vertices)
try {
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  impl->setMeshQuads(meshID, size, vertices);
} catch (...) {
  std::abort();
}

void precicec_setMeshTetrahedron(
    int meshID,
    int firstVertexID,
    int secondVertexID,
    int thirdVertexID,
    int fourthVertexID)
try {
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  impl->setMeshTetrahedron(meshID, firstVertexID, secondVertexID, thirdVertexID, fourthVertexID);
} catch (...) {
  std::abort();
}

void precicec_setMeshTetrahedra(
    int        meshID,
    int        size,
    const int *vertices)
try {
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  impl->setMeshTetrahedra(meshID, size, vertices);
} catch (...) {
  std::abort();
}

void precicec_writeBlockVectorData(
    int           dataID,
    int           size,
    const int *   valueIndices,
    const double *values)
try {
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  impl->writeBlockVectorData(dataID, size, valueIndices, values);
} catch (...) {
  std::abort();
}

void precicec_writeVectorData(
    int           dataID,
    int           valueIndex,
    const double *dataValue)
try {
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  impl->writeVectorData(dataID, valueIndex, dataValue);
} catch (...) {
  std::abort();
}

void precicec_writeBlockScalarData(
    int           dataID,
    int           size,
    const int *   valueIndices,
    const double *values)
try {
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  impl->writeBlockScalarData(dataID, size, valueIndices, values);
} catch (...) {
  std::abort();
}

void precicec_writeScalarData(
    int    dataID,
    int    valueIndex,
    double dataValue)
try {
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  impl->writeScalarData(dataID, valueIndex, dataValue);
} catch (...) {
  std::abort();
}

void precicec_readBlockVectorData(
    int        dataID,
    int        size,
    const int *valueIndices,
    double *   values)
try {
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  impl->readBlockVectorData(dataID, size, valueIndices, values);
} catch (...) {
  std::abort();
}

void precicec_readVectorData(
    int     dataID,
    int     valueIndex,
    double *dataValue)
try {
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  impl->readVectorData(dataID, valueIndex, dataValue);
} catch (...) {
  std::abort();
}

void precicec_readBlockScalarData(
    int        dataID,
    int        size,
    const int *valueIndices,
    double *   values)
try {
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  impl->readBlockScalarData(dataID, size, valueIndices, values);
} catch (...) {
  std::abort();
}

void precicec_readScalarData(
    int     dataID,
    int     valueIndex,
    double *dataValue)
try {
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  impl->readScalarData(dataID, valueIndex, *dataValue);
} catch (...) {
  std::abort();
}

int precicec_requiresGradientDataFor(int dataID)
try {
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  if (impl->requiresGradientDataFor(dataID)) {
    return 1;
  }
  return 0;
} catch (...) {
  std::abort();
}

void precicec_writeScalarGradientData(
    int           dataID,
    int           valueIndex,
    const double *gradientValues)
try {
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  impl->writeScalarGradientData(dataID, valueIndex, gradientValues);
} catch (...) {
  std::abort();
}

void precicec_writeBlockScalarGradientData(
    int           dataID,
    int           size,
    const int *   valueIndices,
    const double *gradientValues)
try {
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  impl->writeBlockScalarGradientData(dataID, size, valueIndices, gradientValues);
} catch (...) {
  std::abort();
}

void precicec_writeVectorGradientData(
    int           dataID,
    int           valueIndex,
    const double *gradientValues)
try {
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  impl->writeVectorGradientData(dataID, valueIndex, gradientValues);
} catch (...) {
  std::abort();
}

void precicec_writeBlockVectorGradientData(
    int           dataID,
    int           size,
    const int *   valueIndices,
    const double *gradientValues)
try {
  PRECICE_CHECK(impl != nullptr, Err, errormsg);
  impl->writeBlockVectorGradientData(dataID, size, valueIndices, gradientValues);
} catch (...) {
  std::abort();
}

const char *precicec_getVersionInformation()
try {
  return precice::versionInformation;
} catch (...) {
  std::abort();
}

void precicec_setMeshAccessRegion(
    const int     meshID,
    const double *boundingBox)
try {
  impl->setMeshAccessRegion(meshID, boundingBox);
} catch (...) {
  std::abort();
}

void precicec_getMeshVerticesAndIDs(
    const int meshID,
    const int size,
    int *     ids,
    double *  coordinates)
try {
  impl->getMeshVerticesAndIDs(meshID, size, ids, coordinates);
} catch (...) {
  std::abort();
}

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif
#ifdef __clang__
#pragma clang diagnostic pop
#endif
