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
  assertion ( impl != nullptr );
  return impl->initialize ();
}

void precicec_initialize_data()
{
  assertion ( impl != nullptr );
  impl->initializeData ();
}

double precicec_advance( double computedTimestepLength )
{
  assertion ( impl != nullptr );
  return impl->advance ( computedTimestepLength );
}

void precicec_finalize()
{
  assertion ( impl != nullptr );
  impl->finalize ();
}

int precicec_getDimensions()
{
  assertion ( impl != nullptr );
  return impl->getDimensions();
}

int precicec_isCouplingOngoing()
{
  assertion ( impl != nullptr );
  if ( impl->isCouplingOngoing() ) {
    return 1;
  }
  return 0;
}

int precicec_isCouplingTimestepComplete()
{
  assertion ( impl != nullptr );
  if ( impl->isTimestepComplete() ){
    return 1;
  }
  return 0;
}

int precicec_isReadDataAvailable()
{
  assertion ( impl != nullptr );
  if ( impl->isReadDataAvailable() ){
     return 1;
  }
  return 0;
}

int precicec_isWriteDataRequired ( double computedTimestepLength )
{
  assertion ( impl != nullptr );
  if ( impl->isWriteDataRequired(computedTimestepLength) ){
     return 1;
  }
  return 0;
}

int precicec_isActionRequired ( const char* action )
{
  assertion ( impl != nullptr );
  assertion ( action != nullptr );
  if ( impl->isActionRequired(std::string(action)) ){
    return 1;
  }
  return 0;
}

void precicec_fulfilledAction ( const char* action )
{
  assertion ( impl != nullptr );
  assertion ( action != nullptr );
  impl->fulfilledAction ( std::string(action) );
}

int precicec_getMeshID ( const char* meshName )
{
  assertion ( impl != nullptr );
  std::string stringMeshName (meshName);
  return impl->getMeshID (stringMeshName);
}

int precicec_hasData ( const char* dataName, int meshID )
{
  assertion ( impl != nullptr );
  std::string stringDataName (dataName);
  return impl->hasData (stringDataName, meshID);
}

int precicec_getDataID ( const char* dataName, int meshID )
{
  assertion ( impl != nullptr );
  std::string stringDataName (dataName);
  return impl->getDataID (stringDataName, meshID);
}

int precicec_setMeshVertex
(
  int           meshID,
  const double* position )
{
  assertion ( impl != nullptr );
  return impl->setMeshVertex ( meshID, position );
}


void precicec_getMeshVertices
(
  int        meshID,
  int        size,
  const int* ids,
  double*    positions )
{
  assertion(impl != nullptr);
  impl->getMeshVertices(meshID, size, ids, positions);
}

void precicec_setMeshVertices
(
  int           meshID,
  int           size,
  const double* positions,
  int*          ids)
{
  assertion(impl != nullptr);
  impl->setMeshVertices(meshID, size, positions, ids );
}

int precicec_getMeshVertexSize
(
  int meshID )
{
  assertion(impl != nullptr);
  return impl->getMeshVertexSize(meshID);
}

int precicec_setMeshEdge
(
  int meshID,
  int firstVertexID,
  int secondVertexID )
{
  assertion ( impl != nullptr );
  return impl->setMeshEdge ( meshID, firstVertexID, secondVertexID );
}

void precicec_setMeshTriangle
(
  int meshID,
  int firstEdgeID,
  int secondEdgeID,
  int thirdEdgeID )
{
  assertion ( impl != nullptr );
  impl->setMeshTriangle ( meshID, firstEdgeID, secondEdgeID, thirdEdgeID );
}

void precicec_setMeshTriangleWithEdges
(
  int meshID,
  int firstVertexID,
  int secondVertexID,
  int thirdVertexID )
{
  assertion ( impl != nullptr );
  impl->setMeshTriangleWithEdges ( meshID, firstVertexID, secondVertexID, thirdVertexID );
}

void precicec_writeBlockVectorData
(
  int           dataID,
  int           size,
  const int*    valueIndices,
  const double* values )
{
  assertion(impl != nullptr);
  impl->writeBlockVectorData(dataID, size, valueIndices, values);
}

void precicec_writeVectorData
(
  int           dataID,
  int           valueIndex,
  const double* dataValue )
{
  assertion ( impl != nullptr );
  impl->writeVectorData ( dataID, valueIndex, dataValue );
}

void precicec_writeBlockScalarData
(
  int           dataID,
  int           size,
  const int*    valueIndices,
  const double* values )
{
  assertion(impl != nullptr);
  impl->writeBlockScalarData(dataID, size, valueIndices, values);
}

void precicec_writeScalarData
(
  int           dataID,
  int           valueIndex,
  const double& dataValue )
{
  assertion ( impl != nullptr );
  impl->writeScalarData ( dataID, valueIndex, dataValue );
}

void precicec_readBlockVectorData
(
  int        dataID,
  int        size,
  const int* valueIndices,
  double*    values )
{
  assertion(impl != nullptr);
  impl->readBlockVectorData(dataID, size, valueIndices, values);
}

void precicec_readVectorData
(
  int     dataID,
  int     valueIndex,
  double* dataValue )
{
  assertion ( impl != nullptr );
  impl->readVectorData (dataID, valueIndex, dataValue);
}

void precicec_readBlockScalarData
(
  int        dataID,
  int        size,
  const int* valueIndices,
  double*    values )
{
  assertion(impl != nullptr);
  impl->readBlockScalarData(dataID, size, valueIndices, values);
}

void precicec_readScalarData
(
  int     dataID,
  int     valueIndex,
  double* dataValue )
{
  assertion ( impl != nullptr );
  impl->readScalarData (dataID, valueIndex, *dataValue);
}

void precicec_mapWriteDataFrom ( int fromMeshID )
{
  assertion ( impl != nullptr );
  impl->mapWriteDataFrom(fromMeshID);
}

void precicec_mapReadDataTo ( int toMeshID )
{
  assertion ( impl != nullptr );
  impl->mapReadDataTo(toMeshID);
}
