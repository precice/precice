extern "C" {
#include "precice/preciceC.h"
}
#include <iostream>
#include <memory>
#include <string>
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

static std::unique_ptr<precice::Participant> impl = nullptr;

static std::string errormsg       = "preCICE has not been created properly. Be sure to call \"precicec_createParticipant\" or \"precicec_createParticipant_withCommunicator\" before any other call to preCICE.";
static std::string errormsgCreate = "preCICE has been created already! Be sure to call either \"precicec_createParticipant\" or \"precicec_createParticipant_withCommunicator\" exactly once.";

void precicec_createParticipant_withCommunicator(
    const char *participantName,
    const char *configFileName,
    int         solverProcessIndex,
    int         solverProcessSize,
    void       *communicator)
try {
  if (impl != nullptr) {
    std::cerr << errormsgCreate;
    std::abort();
  }
  impl = std::make_unique<precice::Participant>(participantName,
                                                configFileName,
                                                solverProcessIndex,
                                                solverProcessSize,
                                                communicator);
} catch (::precice::Error &e) {
  std::abort();
}

void precicec_createParticipant(
    const char *participantName,
    const char *configFileName,
    int         solverProcessIndex,
    int         solverProcessSize)
try {
  if (impl != nullptr) {
    std::cerr << errormsgCreate;
    std::abort();
  }
  impl = std::make_unique<precice::Participant>(participantName,
                                                configFileName,
                                                solverProcessIndex,
                                                solverProcessSize);
} catch (::precice::Error &e) {
  std::abort();
}

void precicec_initialize()
try {
  if (impl == nullptr) {
    std::cerr << errormsg;
    std::abort();
  }
  impl->initialize();
} catch (::precice::Error &e) {
  std::abort();
}

void precicec_advance(double computedTimeStepSize)
try {
  if (impl == nullptr) {
    std::cerr << errormsg;
    std::abort();
  }
  impl->advance(computedTimeStepSize);
} catch (::precice::Error &e) {
  std::abort();
}

void precicec_finalize()
try {
  if (impl == nullptr) {
    std::cerr << errormsg;
    std::abort();
  }
  impl->finalize();
  impl.reset();
} catch (::precice::Error &e) {
  std::abort();
}

int precicec_getMeshDimensions(const char *meshName)
try {
  if (impl == nullptr) {
    std::cerr << errormsg;
    std::abort();
  }
  return impl->getMeshDimensions(meshName);
} catch (::precice::Error &e) {
  std::abort();
  return -1;
}

int precicec_getDataDimensions(const char *meshName, const char *dataName)
try {
  if (impl == nullptr) {
    std::cerr << errormsg;
    std::abort();
  }
  return impl->getDataDimensions(meshName, dataName);
} catch (::precice::Error &e) {
  std::abort();
  return -1;
}

int precicec_isCouplingOngoing()
try {
  if (impl == nullptr) {
    std::cerr << errormsg;
    std::abort();
  }
  if (impl->isCouplingOngoing()) {
    return 1;
  }
  return 0;
} catch (::precice::Error &e) {
  std::abort();
  return -1;
}

int precicec_isTimeWindowComplete()
try {
  if (impl == nullptr) {
    std::cerr << errormsg;
    std::abort();
  }
  if (impl->isTimeWindowComplete()) {
    return 1;
  }
  return 0;
} catch (::precice::Error &e) {
  std::abort();
  return -1;
}

double precicec_getMaxTimeStepSize()
try {
  return impl->getMaxTimeStepSize();
} catch (::precice::Error &e) {
  std::abort();
  return -1.0;
}

int precicec_requiresInitialData()
try {
  if (impl == nullptr) {
    std::cerr << errormsg;
    std::abort();
  }
  return impl->requiresInitialData() ? 1 : 0;
} catch (::precice::Error &e) {
  std::abort();
  return -1;
}

int precicec_requiresWritingCheckpoint()
try {
  if (impl == nullptr) {
    std::cerr << errormsg;
    std::abort();
  }
  return impl->requiresWritingCheckpoint() ? 1 : 0;
} catch (::precice::Error &e) {
  std::abort();
  return -1;
}

int precicec_requiresReadingCheckpoint()
try {
  if (impl == nullptr) {
    std::cerr << errormsg;
    std::abort();
  }
  return impl->requiresReadingCheckpoint() ? 1 : 0;
} catch (::precice::Error &e) {
  std::abort();
  return -1;
}

int precicec_requiresMeshConnectivityFor(const char *meshName)
try {
  if (impl == nullptr) {
    std::cerr << errormsg;
    std::abort();
  }
  if (impl->requiresMeshConnectivityFor(meshName)) {
    return 1;
  }
  return 0;
} catch (::precice::Error &e) {
  std::abort();
  return -1;
}

void precicec_resetMesh(const char *meshName)
try {
  if (impl == nullptr) {
    std::cerr << errormsg;
    std::abort();
  }
  impl->resetMesh(meshName);
} catch (::precice::Error &e) {
  std::abort();
}

int precicec_setMeshVertex(
    const char   *meshName,
    const double *position)
try {
  if (impl == nullptr) {
    std::cerr << errormsg;
    std::abort();
  }
  auto size = static_cast<long unsigned>(impl->getMeshDimensions(meshName));
  return impl->setMeshVertex(meshName, {position, size});
} catch (::precice::Error &e) {
  std::abort();
  return -1;
}

void precicec_setMeshVertices(
    const char   *meshName,
    int           size,
    const double *positions,
    int          *ids)
try {
  if (impl == nullptr) {
    std::cerr << errormsg;
    std::abort();
  }
  auto idsSize = static_cast<long unsigned>(size);
  auto posSize = static_cast<long unsigned>(impl->getMeshDimensions(meshName) * size);
  impl->setMeshVertices(meshName, {positions, posSize}, {ids, idsSize});
} catch (::precice::Error &e) {
  std::abort();
}

int precicec_getMeshVertexSize(
    const char *meshName)
try {
  if (impl == nullptr) {
    std::cerr << errormsg;
    std::abort();
  }
  return impl->getMeshVertexSize(meshName);
} catch (::precice::Error &e) {
  std::abort();
  return -1;
}

void precicec_setMeshEdge(
    const char *meshName,
    int         firstVertexID,
    int         secondVertexID)
try {
  if (impl == nullptr) {
    std::cerr << errormsg;
    std::abort();
  }
  impl->setMeshEdge(meshName, firstVertexID, secondVertexID);
} catch (::precice::Error &e) {
  std::abort();
}

void precicec_setMeshEdges(
    const char *meshName,
    int         size,
    const int  *vertices)
try {
  if (impl == nullptr) {
    std::cerr << errormsg;
    std::abort();
  }
  auto verticesSize = static_cast<long unsigned>(size) * 2;
  impl->setMeshEdges(meshName, {vertices, verticesSize});
} catch (::precice::Error &e) {
  std::abort();
}

void precicec_setMeshTriangle(
    const char *meshName,
    int         firstVertexID,
    int         secondVertexID,
    int         thirdVertexID)
try {
  if (impl == nullptr) {
    std::cerr << errormsg;
    std::abort();
  }
  impl->setMeshTriangle(meshName, firstVertexID, secondVertexID, thirdVertexID);
} catch (::precice::Error &e) {
  std::abort();
}

void precicec_setMeshTriangles(
    const char *meshName,
    int         size,
    const int  *vertices)
try {
  if (impl == nullptr) {
    std::cerr << errormsg;
    std::abort();
  }
  auto verticesSize = static_cast<long unsigned>(size) * 3;
  impl->setMeshTriangles(meshName, {vertices, verticesSize});
} catch (::precice::Error &e) {
  std::abort();
}

void precicec_setMeshQuad(
    const char *meshName,
    int         firstVertexID,
    int         secondVertexID,
    int         thirdVertexID,
    int         fourthVertexID)
try {
  if (impl == nullptr) {
    std::cerr << errormsg;
    std::abort();
  }
  impl->setMeshQuad(meshName, firstVertexID, secondVertexID, thirdVertexID, fourthVertexID);
} catch (::precice::Error &e) {
  std::abort();
}

void precicec_setMeshQuads(
    const char *meshName,
    int         size,
    const int  *vertices)
try {
  if (impl == nullptr) {
    std::cerr << errormsg;
    std::abort();
  }
  auto verticesSize = static_cast<long unsigned>(size) * 4;
  impl->setMeshQuads(meshName, {vertices, verticesSize});
} catch (::precice::Error &e) {
  std::abort();
}

void precicec_setMeshTetrahedron(
    const char *meshName,
    int         firstVertexID,
    int         secondVertexID,
    int         thirdVertexID,
    int         fourthVertexID)
try {
  if (impl == nullptr) {
    std::cerr << errormsg;
    std::abort();
  }
  impl->setMeshTetrahedron(meshName, firstVertexID, secondVertexID, thirdVertexID, fourthVertexID);
} catch (::precice::Error &e) {
  std::abort();
}

void precicec_setMeshTetrahedra(
    const char *meshName,
    int         size,
    const int  *vertices)
try {
  if (impl == nullptr) {
    std::cerr << errormsg;
    std::abort();
  }
  auto verticesSize = static_cast<long unsigned>(size) * 4;
  impl->setMeshTetrahedra(meshName, {vertices, verticesSize});
} catch (::precice::Error &e) {
  std::abort();
}

void precicec_writeData(
    const char   *meshName,
    const char   *dataName,
    int           size,
    const int    *valueIndices,
    const double *values)
try {
  if (impl == nullptr) {
    std::cerr << errormsg;
    std::abort();
  }
  auto dataSize = size * impl->getDataDimensions(meshName, dataName);
  impl->writeData(meshName, dataName, {valueIndices, static_cast<unsigned long>(size)}, {values, static_cast<unsigned long>(dataSize)});
} catch (::precice::Error &e) {
  std::abort();
}

void precicec_readData(
    const char *meshName,
    const char *dataName,
    int         size,
    const int  *valueIndices,
    double      relativeReadTime,
    double     *values)
try {
  if (impl == nullptr) {
    std::cerr << errormsg;
    std::abort();
  }
  auto dataSize = size * impl->getDataDimensions(meshName, dataName);
  impl->readData(meshName, dataName, {valueIndices, static_cast<unsigned long>(size)}, relativeReadTime, {values, static_cast<unsigned long>(dataSize)});
} catch (::precice::Error &e) {
  std::abort();
}

void precicec_writeAndMapData(
    const char   *meshName,
    const char   *dataName,
    int           size,
    const double *coordinates,
    const double *values)
try {
  if (impl == nullptr) {
    std::cerr << errormsg;
    std::abort();
  }
  auto coordinatesSize = size * impl->getMeshDimensions(meshName);
  auto dataSize        = size * impl->getDataDimensions(meshName, dataName);
  impl->writeAndMapData(meshName, dataName, {coordinates, static_cast<unsigned long>(coordinatesSize)}, {values, static_cast<unsigned long>(dataSize)});
} catch (::precice::Error &e) {
  std::abort();
}

void precicec_mapAndReadData(
    const char   *meshName,
    const char   *dataName,
    int           size,
    const double *coordinates,
    double        relativeReadTime,
    double       *values)
try {
  if (impl == nullptr) {
    std::cerr << errormsg;
    std::abort();
  }
  auto coordinatesSize = size * impl->getMeshDimensions(meshName);
  auto dataSize        = size * impl->getDataDimensions(meshName, dataName);
  impl->mapAndReadData(meshName, dataName, {coordinates, static_cast<unsigned long>(coordinatesSize)}, relativeReadTime, {values, static_cast<unsigned long>(dataSize)});
} catch (::precice::Error &e) {
  std::abort();
}

int precicec_requiresGradientDataFor(const char *meshName,
                                     const char *dataName)
try {
  if (impl == nullptr) {
    std::cerr << errormsg;
    std::abort();
  }
  if (impl->requiresGradientDataFor(meshName, dataName)) {
    return 1;
  }
  return 0;
} catch (::precice::Error &e) {
  std::abort();
  return -1;
}

void precicec_writeGradientData(
    const char   *meshName,
    const char   *dataName,
    int           size,
    const int    *valueIndices,
    const double *gradients)
try {
  if (impl == nullptr) {
    std::cerr << errormsg;
    std::abort();
  }
  auto gradientComponents = impl->getDataDimensions(meshName, dataName) * impl->getMeshDimensions(meshName);
  auto gradientSize       = size * gradientComponents;
  impl->writeGradientData(meshName, dataName, {valueIndices, static_cast<unsigned long>(size)}, {gradients, static_cast<unsigned long>(gradientSize)});
} catch (::precice::Error &e) {
  std::abort();
}

const char *precicec_getVersionInformation()
{
  return precice::versionInformation;
}

void precicec_setMeshAccessRegion(
    const char   *meshName,
    const double *boundingBox)
try {
  auto bbSize = static_cast<long unsigned>(impl->getMeshDimensions(meshName)) * 2;
  impl->setMeshAccessRegion(meshName, {boundingBox, bbSize});
} catch (::precice::Error &e) {
  std::abort();
}

void precicec_getMeshVertexIDsAndCoordinates(
    const char *meshName,
    const int   size,
    int        *ids,
    double     *coordinates)
try {
  auto coordinatesSize = static_cast<long unsigned>(impl->getMeshDimensions(meshName) * size);
  impl->getMeshVertexIDsAndCoordinates(meshName, {ids, static_cast<unsigned long>(size)}, {coordinates, coordinatesSize});
} catch (::precice::Error &e) {
  std::abort();
}

void precicec_startProfilingSection(const char *sectionName)
try {
  impl->startProfilingSection(sectionName);
} catch (::precice::Error &e) {
  std::abort();
}

void precicec_stopLastProfilingSection()
try {
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
