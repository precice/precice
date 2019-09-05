extern "C" {
#include "SolverInterfaceC.h"
}
#include "precice/impl/SolverInterfaceImpl.hpp"
#include "utils/assertion.hpp"
#include <string>

static precice::impl::SolverInterfaceImpl* impl = nullptr;

void precicec_createSolverInterface
(
  const char* participantName,
  const char* configFileName,
  int         solverProcessIndex,
  int         solverProcessSize )
{
  std::string stringAccessorName ( participantName );
  impl = new precice::impl::SolverInterfaceImpl ( stringAccessorName,
      solverProcessIndex, solverProcessSize, false );
  std::string stringConfigFileName ( configFileName );
  impl->configure ( stringConfigFileName );
}

double precicec_initialize()
{
  PRECICE_ASSERT( impl != nullptr );
  return impl->initialize ();
}

void precicec_initialize_data()
{
  PRECICE_ASSERT( impl != nullptr );
  impl->initializeData ();
}

double precicec_advance( double computedTimestepLength )
{
  PRECICE_ASSERT( impl != nullptr );
  return impl->advance ( computedTimestepLength );
}

void precicec_finalize()
{
  PRECICE_ASSERT( impl != nullptr );
  impl->finalize ();
}

int precicec_getDimensions()
{
  PRECICE_ASSERT( impl != nullptr );
  return impl->getDimensions();
}

int precicec_isCouplingOngoing()
{
  PRECICE_ASSERT( impl != nullptr );
  if ( impl->isCouplingOngoing() ) {
    return 1;
  }
  return 0;
}

int precicec_isCouplingTimestepComplete()
{
  PRECICE_ASSERT( impl != nullptr );
  if ( impl->isTimestepComplete() ){
    return 1;
  }
  return 0;
}

int precicec_hasToEvaluateSurrogateModel()
{
  PRECICE_ASSERT( impl != nullptr);
  if (impl->hasToEvaluateSurrogateModel() ){
    return 1;
  }
  return 0;
}

int precicec_hasToEvaluateFineModel()
{
  PRECICE_ASSERT( impl != nullptr);
  if (impl->hasToEvaluateFineModel() ){
    return 1;
  }
  return 0;
}

int precicec_isReadDataAvailable()
{
  PRECICE_ASSERT( impl != nullptr );
  if ( impl->isReadDataAvailable() ){
     return 1;
  }
  return 0;
}

int precicec_isWriteDataRequired ( double computedTimestepLength )
{
  PRECICE_ASSERT( impl != nullptr );
  if ( impl->isWriteDataRequired(computedTimestepLength) ){
     return 1;
  }
  return 0;
}

int precicec_isActionRequired ( const char* action )
{
  PRECICE_ASSERT( impl != nullptr );
  PRECICE_ASSERT( action != nullptr );
  if ( impl->isActionRequired(std::string(action)) ){
    return 1;
  }
  return 0;
}

void precicec_fulfilledAction ( const char* action )
{
  PRECICE_ASSERT( impl != nullptr );
  PRECICE_ASSERT( action != nullptr );
  impl->fulfilledAction ( std::string(action) );
}

int precicec_hasMesh ( const char* meshName)
{
  PRECICE_ASSERT( impl != nullptr );
  std::string stringMeshName (meshName);
  if ( impl->hasMesh (stringMeshName) ){
    return 1;
  }
  return 0;
}

int precicec_getMeshID ( const char* meshName )
{
  PRECICE_ASSERT( impl != nullptr );
  std::string stringMeshName (meshName);
  return impl->getMeshID (stringMeshName);
}

int precicec_hasData ( const char* dataName, int meshID )
{
  PRECICE_ASSERT( impl != nullptr );
  std::string stringDataName (dataName);
  return impl->hasData (stringDataName, meshID);
}

int precicec_getDataID ( const char* dataName, int meshID )
{
  PRECICE_ASSERT( impl != nullptr );
  std::string stringDataName (dataName);
  return impl->getDataID (stringDataName, meshID);
}

int precicec_setMeshVertex
(
  int           meshID,
  const double* position )
{
  PRECICE_ASSERT( impl != nullptr );
  return impl->setMeshVertex ( meshID, position );
}


void precicec_getMeshVertices
(
  int        meshID,
  int        size,
  const int* ids,
  double*    positions )
{
  PRECICE_ASSERT(impl != nullptr);
  impl->getMeshVertices(meshID, size, ids, positions);
}

void precicec_setMeshVertices
(
  int           meshID,
  int           size,
  const double* positions,
  int*          ids)
{
  PRECICE_ASSERT(impl != nullptr);
  impl->setMeshVertices(meshID, size, positions, ids );
}

int precicec_getMeshVertexSize
(
  int meshID )
{
  PRECICE_ASSERT(impl != nullptr);
  return impl->getMeshVertexSize(meshID);
}

void precicec_getMeshVertexIDsFromPositions
(
  int           meshID,
  int           size,
  const double* positions,
  int*          ids)
{
  PRECICE_ASSERT(impl != nullptr);
  impl->getMeshVertexIDsFromPositions(meshID,size,positions,ids);
}

int precicec_setMeshEdge
(
  int meshID,
  int firstVertexID,
  int secondVertexID )
{
  PRECICE_ASSERT( impl != nullptr );
  return impl->setMeshEdge ( meshID, firstVertexID, secondVertexID );
}

void precicec_setMeshTriangle
(
  int meshID,
  int firstEdgeID,
  int secondEdgeID,
  int thirdEdgeID )
{
  PRECICE_ASSERT( impl != nullptr );
  impl->setMeshTriangle ( meshID, firstEdgeID, secondEdgeID, thirdEdgeID );
}

void precicec_setMeshTriangleWithEdges
(
  int meshID,
  int firstVertexID,
  int secondVertexID,
  int thirdVertexID )
{
  PRECICE_ASSERT( impl != nullptr );
  impl->setMeshTriangleWithEdges ( meshID, firstVertexID, secondVertexID, thirdVertexID );
}

void precicec_setMeshQuad
(
  int meshID,
  int firstEdgeID,
  int secondEdgeID,
  int thirdEdgeID,
  int fourthEdgeID )
{
  PRECICE_ASSERT( impl != nullptr );
  impl->setMeshQuad(meshID,firstEdgeID,secondEdgeID,thirdEdgeID,fourthEdgeID);
}

void precicec_setMeshQuadWithEdges
(
  int meshID,
  int firstVertexID,
  int secondVertexID,
  int thirdVertexID,
  int fourthVertexID )
{
  PRECICE_ASSERT( impl != nullptr );
  impl->setMeshQuadWithEdges(meshID,firstVertexID,secondVertexID,thirdVertexID,fourthVertexID);
}

void precicec_writeBlockVectorData
(
  int           dataID,
  int           size,
  const int*    valueIndices,
  const double* values )
{
  PRECICE_ASSERT(impl != nullptr);
  impl->writeBlockVectorData(dataID, size, valueIndices, values);
}

void precicec_writeVectorData
(
  int           dataID,
  int           valueIndex,
  const double* dataValue )
{
  PRECICE_ASSERT( impl != nullptr );
  impl->writeVectorData ( dataID, valueIndex, dataValue );
}

void precicec_writeBlockScalarData
(
  int           dataID,
  int           size,
  const int*    valueIndices,
  const double* values )
{
  PRECICE_ASSERT(impl != nullptr);
  impl->writeBlockScalarData(dataID, size, valueIndices, values);
}

void precicec_writeScalarData
(
  int    dataID,
  int    valueIndex,
  double dataValue )
{
  PRECICE_ASSERT( impl != nullptr );
  impl->writeScalarData ( dataID, valueIndex, dataValue );
}

void precicec_readBlockVectorData
(
  int        dataID,
  int        size,
  const int* valueIndices,
  double*    values )
{
  PRECICE_ASSERT(impl != nullptr);
  impl->readBlockVectorData(dataID, size, valueIndices, values);
}

void precicec_readVectorData
(
  int     dataID,
  int     valueIndex,
  double* dataValue )
{
  PRECICE_ASSERT( impl != nullptr );
  impl->readVectorData (dataID, valueIndex, dataValue);
}

void precicec_readBlockScalarData
(
  int        dataID,
  int        size,
  const int* valueIndices,
  double*    values )
{
  PRECICE_ASSERT(impl != nullptr);
  impl->readBlockScalarData(dataID, size, valueIndices, values);
}

void precicec_readScalarData
(
  int     dataID,
  int     valueIndex,
  double* dataValue )
{
  PRECICE_ASSERT( impl != nullptr );
  impl->readScalarData (dataID, valueIndex, *dataValue);
}

void precicec_mapWriteDataFrom ( int fromMeshID )
{
  PRECICE_ASSERT( impl != nullptr );
  impl->mapWriteDataFrom(fromMeshID);
}

void precicec_mapReadDataTo ( int toMeshID )
{
  PRECICE_ASSERT( impl != nullptr );
  impl->mapReadDataTo(toMeshID);
}

const char* precicec_actionWriteInitialData()
{
  return precice::constants::actionWriteInitialData().c_str();
}

const char* precicec_actionWriteIterationCheckpoint()
{
  return precice::constants::actionWriteIterationCheckpoint().c_str();
}

const char* precicec_actionReadIterationCheckpoint()
{
  return precice::constants::actionReadIterationCheckpoint().c_str();
}
