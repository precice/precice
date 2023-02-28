#include "precice/SolverInterface.hpp"
#include "cplscheme/Constants.hpp"
#include "precice/impl/SolverInterfaceImpl.hpp"
#include "precice/impl/versions.hpp"

namespace precice {

SolverInterface::SolverInterface(
    const std::string &participantName,
    const std::string &configurationFileName,
    int                solverProcessIndex,
    int                solverProcessSize)
    : _impl(new impl::SolverInterfaceImpl(participantName, configurationFileName, solverProcessIndex, solverProcessSize))
{
}

SolverInterface::SolverInterface(
    const std::string &participantName,
    const std::string &configurationFileName,
    int                solverProcessIndex,
    int                solverProcessSize,
    void *             communicator)
    : _impl(new impl::SolverInterfaceImpl(participantName, configurationFileName, solverProcessIndex, solverProcessSize, communicator))
{
}

SolverInterface::~SolverInterface() = default;

double SolverInterface::initialize()
{
  return _impl->initialize();
}

double SolverInterface::advance(
    double computedTimestepLength)
{
  return _impl->advance(computedTimestepLength);
}

void SolverInterface::finalize()
{
  return _impl->finalize();
}

int SolverInterface::getDimensions() const
{
  return _impl->getDimensions();
}

bool SolverInterface::isCouplingOngoing() const
{
  return _impl->isCouplingOngoing();
}

bool SolverInterface::isTimeWindowComplete() const
{
  return _impl->isTimeWindowComplete();
}

bool SolverInterface::requiresInitialData()
{
  return _impl->requiresInitialData();
}

bool SolverInterface::requiresReadingCheckpoint()
{
  return _impl->requiresReadingCheckpoint();
}

bool SolverInterface::requiresWritingCheckpoint()
{
  return _impl->requiresWritingCheckpoint();
}

bool SolverInterface::hasMesh(
    const std::string &meshName) const
{
  return _impl->hasMesh(meshName);
}

int SolverInterface::getMeshID(
    const std::string &meshName) const
{
  return _impl->getMeshID(meshName);
}

bool SolverInterface::requiresMeshConnectivityFor(int meshID) const
{
  return _impl->requiresMeshConnectivityFor(meshID);
}

bool SolverInterface::requiresGradientDataFor(int dataID) const
{
  return _impl->requiresGradientDataFor(dataID);
}

bool SolverInterface::hasData(
    const std::string &dataName, int meshID) const
{
  return _impl->hasData(dataName, meshID);
}

int SolverInterface::getDataID(
    const std::string &dataName, int meshID) const
{
  return _impl->getDataID(dataName, meshID);
}

//void SolverInterface:: resetMesh
//(
//  int meshID )
//{
//  _impl->resetMesh(meshID);
//}

int SolverInterface::setMeshVertex(
    int           meshID,
    const double *position)
{
  return _impl->setMeshVertex(meshID, position);
}

int SolverInterface::getMeshVertexSize(
    int meshID) const
{
  return _impl->getMeshVertexSize(meshID);
}

void SolverInterface::setMeshVertices(
    int           meshID,
    int           size,
    const double *positions,
    int *         ids)
{
  _impl->setMeshVertices(meshID, size, positions, ids);
}

void SolverInterface::setMeshEdge(
    int meshID,
    int firstVertexID,
    int secondVertexID)
{
  _impl->setMeshEdge(meshID, firstVertexID, secondVertexID);
}

void SolverInterface::setMeshEdges(
    int        meshID,
    int        size,
    const int *vertices)
{
  _impl->setMeshEdges(meshID, size, vertices);
}

void SolverInterface::setMeshTriangle(
    int meshID,
    int firstVertexID,
    int secondVertexID,
    int thirdVertexID)
{
  _impl->setMeshTriangle(meshID, firstVertexID, secondVertexID, thirdVertexID);
}

void SolverInterface::setMeshTriangles(
    int        meshID,
    int        size,
    const int *vertices)
{
  _impl->setMeshTriangles(meshID, size, vertices);
}

void SolverInterface::setMeshQuad(
    int meshID,
    int firstVertexID,
    int secondVertexID,
    int thirdVertexID,
    int fourthVertexID)
{
  _impl->setMeshQuad(meshID, firstVertexID, secondVertexID, thirdVertexID,
                     fourthVertexID);
}

void SolverInterface::setMeshQuads(
    int        meshID,
    int        size,
    const int *vertices)
{
  _impl->setMeshQuads(meshID, size, vertices);
}

void SolverInterface::setMeshTetrahedron(
    int meshID,
    int firstVertexID,
    int secondVertexID,
    int thirdVertexID,
    int fourthVertexID)
{
  _impl->setMeshTetrahedron(meshID, firstVertexID, secondVertexID, thirdVertexID,
                            fourthVertexID);
}

void SolverInterface::setMeshTetrahedra(
    int        meshID,
    int        size,
    const int *vertices)
{
  _impl->setMeshTetrahedra(meshID, size, vertices);
}

void SolverInterface::writeBlockVectorData(
    int           dataID,
    int           size,
    const int *   valueIndices,
    const double *values)
{
  _impl->writeBlockVectorData(dataID, size, valueIndices, values);
}

void SolverInterface::writeBlockVectorGradientData(
    int           dataID,
    int           size,
    const int *   valueIndices,
    const double *gradientValues)
{
  _impl->writeBlockVectorGradientData(dataID, size, valueIndices, gradientValues);
}

void SolverInterface::writeVectorData(
    int           dataID,
    int           valueIndex,
    const double *value)
{
  _impl->writeVectorData(dataID, valueIndex, value);
}

void SolverInterface::writeVectorGradientData(
    int           dataID,
    int           valueIndex,
    const double *gradientValues)
{
  _impl->writeVectorGradientData(dataID, valueIndex, gradientValues);
}

void SolverInterface::writeBlockScalarData(
    int           dataID,
    int           size,
    const int *   valueIndices,
    const double *values)
{
  _impl->writeBlockScalarData(dataID, size, valueIndices, values);
}

void SolverInterface::writeBlockScalarGradientData(
    int           dataID,
    int           size,
    const int *   valueIndices,
    const double *gradientValues)
{
  _impl->writeBlockScalarGradientData(dataID, size, valueIndices, gradientValues);
}

void SolverInterface::writeScalarData(
    int    dataID,
    int    valueIndex,
    double value)
{
  _impl->writeScalarData(dataID, valueIndex, value);
}

void SolverInterface::writeScalarGradientData(
    int           dataID,
    int           valueIndex,
    const double *gradientValues)
{
  _impl->writeScalarGradientData(dataID, valueIndex, gradientValues);
}

void SolverInterface::readBlockVectorData(
    int        dataID,
    int        size,
    const int *valueIndices,
    double *   values) const
{
  _impl->readBlockVectorData(dataID, size, valueIndices, values);
}

void SolverInterface::readBlockVectorData(
    int        dataID,
    int        size,
    const int *valueIndices,
    double     relativeReadTime,
    double *   values) const
{
  _impl->readBlockVectorData(dataID, size, valueIndices, relativeReadTime, values);
}

void SolverInterface::readVectorData(
    int     dataID,
    int     valueIndex,
    double *value) const
{
  _impl->readVectorData(dataID, valueIndex, value);
}

void SolverInterface::readVectorData(
    int     dataID,
    int     valueIndex,
    double  relativeReadTime,
    double *value) const
{
  // @todo: needs testing!
  _impl->readVectorData(dataID, valueIndex, relativeReadTime, value);
}

void SolverInterface::readBlockScalarData(
    int        dataID,
    int        size,
    const int *valueIndices,
    double *   values) const
{
  _impl->readBlockScalarData(dataID, size, valueIndices, values);
}

void SolverInterface::readBlockScalarData(
    int        dataID,
    int        size,
    const int *valueIndices,
    double     relativeReadTime,
    double *   values) const
{
  _impl->readBlockScalarData(dataID, size, valueIndices, relativeReadTime, values);
}

void SolverInterface::readScalarData(
    int     dataID,
    int     valueIndex,
    double &value) const
{
  _impl->readScalarData(dataID, valueIndex, value);
}

void SolverInterface::readScalarData(
    int     dataID,
    int     valueIndex,
    double  relativeReadTime,
    double &value) const
{
  _impl->readScalarData(dataID, valueIndex, relativeReadTime, value);
}

void SolverInterface::setMeshAccessRegion(const int     meshID,
                                          const double *boundingBox) const
{
  _impl->setMeshAccessRegion(meshID, boundingBox);
}

void SolverInterface::getMeshVerticesAndIDs(const int meshID,
                                            const int size,
                                            int *     ids,
                                            double *  coordinates) const
{
  _impl->getMeshVerticesAndIDs(meshID, size, ids, coordinates);
}

} // namespace precice
