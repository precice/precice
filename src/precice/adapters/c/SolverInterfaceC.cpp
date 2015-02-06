// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
extern "C" {
#include "SolverInterfaceC.h"
}
#include "precice/impl/SolverInterfaceImpl.hpp"
#include "utils/Globals.hpp"
#include "utils/Dimensions.hpp"
#include <string>
#include "boost/smart_ptr.hpp"
#include <vector>

static precice::impl::SolverInterfaceImpl* impl = NULL;

void precicec_createSolverInterface
(
  const char* accessorName,
  const char* configFileName,
  int         solverProcessIndex,
  int         solverProcessSize )
{
  std::string stringAccessorName ( accessorName );
  impl = new precice::impl::SolverInterfaceImpl ( stringAccessorName,
      solverProcessIndex, solverProcessSize, false );
  std::string stringConfigFileName ( configFileName );
  impl->configure ( stringConfigFileName );
}

double precicec_initialize()
{
  assertion ( impl != NULL );
  return impl->initialize ();
}

double precicec_advance( double computedTimestepLength )
{
  assertion ( impl != NULL );
  return impl->advance ( computedTimestepLength );
}

void precicec_finalize()
{
  assertion ( impl != NULL );
  impl->finalize ();
}

int precicec_getDimensions()
{
  assertion ( impl != NULL );
  return impl->getDimensions();
}

int precicec_isCouplingOngoing()
{
  assertion ( impl != NULL );
  if ( impl->isCouplingOngoing() ) {
    return 1;
  }
  return 0;
}

int precicec_isCouplingTimestepComplete()
{
  assertion ( impl != NULL );
  if ( impl->isTimestepComplete() ){
    return 1;
  }
  return 0;
}

int precicec_isReadDataAvailable()
{
  assertion ( impl != NULL );
  if ( impl->isReadDataAvailable() ){
     return 1;
  }
  return 0;
}

int precicec_isWriteDataRequired ( double computedTimestepLength )
{
  assertion ( impl != NULL );
  if ( impl->isWriteDataRequired(computedTimestepLength) ){
     return 1;
  }
  return 0;
}

int precicec_isActionRequired ( const char* action )
{
  assertion ( impl != NULL );
  assertion ( action != NULL );
  if ( impl->isActionRequired(std::string(action)) ){
    return 1;
  }
  return 0;
}

void precicec_fulfilledAction ( const char* action )
{
  assertion ( impl != NULL );
  assertion ( action != NULL );
  impl->fulfilledAction ( std::string(action) );
}

int precicec_getMeshID ( const char* geometryName )
{
  assertion ( impl != NULL );
  std::string stringGeometryName (geometryName);
  return impl->getMeshID (stringGeometryName);
}

int precicec_hasData ( const char* dataName, int meshID )
{
  assertion ( impl != NULL );
  std::string stringDataName (dataName);
  return impl->hasData (stringDataName, meshID);
}

int precicec_getDataID ( const char* dataName, int meshID )
{
  assertion ( impl != NULL );
  std::string stringDataName (dataName);
  return impl->getDataID (stringDataName, meshID);
}

int precicec_setMeshVertex
(
  int           meshID,
  const double* position )
{
  assertion ( impl != NULL );
  return impl->setMeshVertex ( meshID, position );
}


void precicec_getMeshVertices
(
  int     meshID,
  int     size,
  int*    ids,
  double* positions )
{
  assertion(impl != NULL);
  impl->getMeshVertices(meshID, size, ids, positions);
}


int precicec_getMeshVertexSize
(
  int meshID )
{
  assertion(impl != NULL);
  return impl->getMeshVertexSize(meshID);
}

int precicec_setMeshEdge
(
  int meshID,
  int firstVertexID,
  int secondVertexID )
{
  assertion ( impl != NULL );
  return impl->setMeshEdge ( meshID, firstVertexID, secondVertexID );
}

void precicec_setMeshTriangle
(
  int meshID,
  int firstEdgeID,
  int secondEdgeID,
  int thirdEdgeID )
{
  assertion ( impl != NULL );
  impl->setMeshTriangle ( meshID, firstEdgeID, secondEdgeID, thirdEdgeID );
}

void precicec_setMeshTriangleWithEdges
(
  int meshID,
  int firstVertexID,
  int secondVertexID,
  int thirdVertexID )
{
  assertion ( impl != NULL );
  impl->setMeshTriangleWithEdges ( meshID, firstVertexID, secondVertexID, thirdVertexID );
}

void precicec_writeBlockVectorData
(
  int     dataID,
  int     size,
  int*    valueIndices,
  double* values )
{
  assertion(impl != NULL);
  impl->writeBlockVectorData(dataID, size, valueIndices, values);
}

void precicec_writeVectorData
(
  int           dataID,
  int           valueIndex,
  const double* dataValue )
{
  assertion ( impl != NULL );
  impl->writeVectorData ( dataID, valueIndex, dataValue );
}

void precicec_writeBlockScalarData
(
  int     dataID,
  int     size,
  int*    valueIndices,
  double* values )
{
  assertion(impl != NULL);
  impl->writeBlockScalarData(dataID, size, valueIndices, values);
}

void precicec_writeScalarData
(
  int    dataID,
  int    valueIndex,
  double dataValue )
{
  assertion ( impl != NULL );
  impl->writeScalarData ( dataID, valueIndex, dataValue );
}

void precicec_readBlockVectorData
(
  int     dataID,
  int     size,
  int*    valueIndices,
  double* values )
{
  assertion(impl != NULL);
  impl->readBlockVectorData(dataID, size, valueIndices, values);
}

void precicec_readVectorData
(
  int     dataID,
  int     valueIndex,
  double* dataValue )
{
  assertion ( impl != NULL );
  impl->readVectorData (dataID, valueIndex, dataValue);
}

void precicec_readBlockScalarData
(
  int     dataID,
  int     size,
  int*    valueIndices,
  double* values )
{
  assertion(impl != NULL);
  impl->readBlockScalarData(dataID, size, valueIndices, values);
}

void precicec_readScalarData
(
  int     dataID,
  int     valueIndex,
  double* dataValue )
{
  assertion ( impl != NULL );
  impl->readScalarData (dataID, valueIndex, *dataValue);
}

void precicec_mapWriteDataFrom ( int fromMeshID )
{
  assertion ( impl != NULL );
  impl->mapWriteDataFrom(fromMeshID);
}

void precicec_mapReadDataTo ( int toMeshID )
{
  assertion ( impl != NULL );
  impl->mapReadDataTo(toMeshID);
}

void precicec_exportMesh
(
  const char* filenameSuffix )
{
  assertion ( impl != NULL );
  std::string stringFilenameSuffix ( filenameSuffix );
  impl->exportMesh ( stringFilenameSuffix );
}
