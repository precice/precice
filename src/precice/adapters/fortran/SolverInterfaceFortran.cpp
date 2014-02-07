#include "SolverInterfaceFortran.hpp"
#include "precice/impl/SolverInterfaceImpl.hpp"
#include <iostream>
#include <string>

using namespace std;

static precice::impl::SolverInterfaceImpl* impl = NULL;

namespace precice {
  namespace impl {
    /**
     * @brief Returns length of string without trailing whitespaces.
     */
    int strippedLength( const char* string, int length );
  }
}

void precicef_create_
(
  const char* accessorName,
  const char* configFileName,
  const int*  solverProcessIndex,
  const int*  solverProcessSize,
  int   lengthAccessorName,
  int   lengthConfigFileName )
{
  //cout << "lengthAccessorName: " << lengthAccessorName << endl;
  //cout << "lengthConfigFileName: " << lengthConfigFileName << endl;
  //cout << "solverProcessIndex: " << *solverProcessIndex << endl;
  //cout << "solverProcessSize: " << *solverProcessSize << endl;
  int strippedLength = precice::impl::strippedLength(accessorName,lengthAccessorName);
  string stringAccessorName(accessorName, strippedLength);
  strippedLength = precice::impl::strippedLength(configFileName,lengthConfigFileName);
  string stringConfigFileName(configFileName, strippedLength);
  //cout << "Accessor: " << stringAccessorName << "!" << endl;
  //cout << "Config  : " << stringConfigFileName << "!" << endl;
  impl = new precice::impl::SolverInterfaceImpl (stringAccessorName,
         *solverProcessIndex, *solverProcessSize, false);
  impl->configure(stringConfigFileName);
}

void precicef_initialize_
(
  double* timestepLengthLimit)
{
  assertion(impl != NULL);
  *timestepLengthLimit = impl->initialize();
}

void precicef_initialize_data_()
{
  assertion(impl != NULL);
  impl->initializeData();
}

void precicef_advance_
(
  double* timestepLengthLimit)
{
  assertion(impl != NULL);
  *timestepLengthLimit = impl->advance(*timestepLengthLimit);
}

void precicef_finalize_()
{
  assertion(impl != NULL);
  impl->finalize();
}

void precicef_get_dims_
(
  int* dimensions )
{
  assertion(impl != NULL);
  *dimensions = impl->getDimensions();
}

void precicef_ongoing_
(
  int* isOngoing )
{
  assertion(impl != NULL);
  if (impl->isCouplingOngoing()){
    *isOngoing = 1;
  }
  else {
    *isOngoing = 0;
  }
}

void precicef_write_data_required_
(
  double* computedTimestepLength,
  int*    isRequired )
{
  assertion(impl != NULL);
  if (impl->isWriteDataRequired(*computedTimestepLength)){
    *isRequired = 1;
  }
  else {
    *isRequired = 0;
  }
}

void precicef_read_data_available_
(
  int* isAvailable )
{
  assertion(impl != NULL);
  if (impl->isReadDataAvailable()){
    *isAvailable = 1;
  }
  else {
    *isAvailable = 0;
  }
}

void precicef_action_required_
(
  const char* action,
  int*        isRequired,
  int         lengthAction )
{
  assertion(impl != NULL);
  //assertion(lengthAction > 1);
  //std::cout << "lengthAction: " << lengthAction << std::endl;
  //std::cout << "Action:";
  //for (int i=0; i < lengthAction; i++){
  //  std::cout << " a[" << i << "]=\"" << action[i] << "\"";
  //}
  //std::cout << std::endl;
  int strippedLength = precice::impl::strippedLength(action, lengthAction);
  //std::cout << "strippedLength: " << strippedLength << std::endl;
  //assertion(strippedLength > 1);
  string stringAction(action, strippedLength);
  if (impl->isActionRequired(stringAction)){
    *isRequired = 1;
  }
  else {
    *isRequired = 0;
  }
}

void precicef_fulfilled_action_
(
  const char* action,
  int         lengthAction )
{
  assertion(impl != NULL);
  int strippedLength = precice::impl::strippedLength(action, lengthAction);
  string stringAction(action, strippedLength);
  impl->fulfilledAction(stringAction);
}

void precicef_get_mesh_id_
(
  const char* geometryName,
  int*        meshID,
  int         lengthGeometryName )
{
  assertion(impl != NULL);
  int strippedLength = precice::impl::strippedLength(geometryName, lengthGeometryName);
  string stringGeometryName(geometryName, strippedLength);
  *meshID = impl->getMeshID(stringGeometryName);
}

void precicef_has_data_
(
  const char* dataName,
  int*        hasData,
  int         lengthDataName )
{
  assertion(impl != NULL);
  int strippedLength = precice::impl::strippedLength(dataName, lengthDataName);
  string stringDataName(dataName, strippedLength);
  if (impl->hasData(stringDataName)){
    *hasData = 1;
  }
  else {
    *hasData = 0;
  }
}

void precicef_get_data_id_
(
  const char* dataName,
  int*        dataID,
  int         lengthDataName )
{
  assertion(impl != NULL);
  int strippedLength = precice::impl::strippedLength(dataName, lengthDataName);
  string stringDataName(dataName, strippedLength);
  *dataID = impl->getDataID(stringDataName);
}

void precicef_set_vertex_
(
  const int*    meshID,
  const double* position,
  int*          vertexID )
{
  assertion(impl != NULL);
  *vertexID = impl->setMeshVertex(*meshID, position);
}

void precicef_set_read_pos_
(
  const int*    meshID,
  const double* position,
  int*          positionID )
{
  assertion(impl != NULL);
  *positionID = impl->setReadPosition(*meshID, position);
}

void precicef_set_read_poss_
(
  const int*    meshID,
  const int*    size,
  double*       positions,
  int*          positionIDs )
{
  assertion(impl != NULL);
  impl->setReadPositions(*meshID, *size, positions, positionIDs);
}

void precicef_set_write_pos_
(
  const int*    meshID,
  const double* position,
  int*          positionID )
{
  assertion(impl != NULL);
  *positionID = impl->setWritePosition(*meshID, position);
}

void precicef_set_write_poss_
(
  const int*    meshID,
  const int*    size,
  double*       positions,
  int*          positionIDs )
{
  assertion(impl != NULL);
  impl->setWritePositions(*meshID, *size, positions, positionIDs);
}

void precicef_set_edge_
(
  const int* meshID,
  const int* firstVertexID,
  const int* secondVertexID,
  int* edgeID )
{
  assertion(impl != NULL);
  *edgeID = impl->setMeshEdge(*meshID, *firstVertexID, *secondVertexID);
}

void precicef_set_triangle_
(
  const int* meshID,
  const int* firstEdgeID,
  const int* secondEdgeID,
  const int* thirdEdgeID )
{
  assertion(impl != NULL);
  impl->setMeshTriangle(*meshID, *firstEdgeID, *secondEdgeID, *thirdEdgeID);
}

void precicef_set_triangle_we_
(
  const int* meshID,
  const int* firstVertexID,
  const int* secondVertexID,
  const int* thirdVertexID )
{
  assertion(impl != NULL);
  impl->setMeshTriangleWithEdges(*meshID, *firstVertexID, *secondVertexID, *thirdVertexID);
}

void precicef_write_bvdata_
(
  const int* dataID,
  const int* size,
  int*       valueIndices,
  double*    values )
{
  assertion(impl != NULL);
  impl->writeBlockVectorData(*dataID, *size, valueIndices, values);
}

void precicef_write_vdata_
(
  const int*    dataID,
  const int*    valueIndex,
  const double* dataValue )
{
  assertion(impl != NULL);
  impl->writeVectorData(*dataID, *valueIndex, dataValue);
}

void precicef_write_bsdata_
(
  const int* dataID,
  const int* size,
  int*       valueIndices,
  double*    values )
{
  assertion(impl != NULL);
  impl->writeBlockScalarData(*dataID, *size, valueIndices, values);
}

void precicef_write_sdata_
(
  const int*    dataID,
  const int*    valueIndex,
  const double* dataValue )
{
  assertion(impl != NULL);
  impl->writeScalarData(*dataID, *valueIndex, *dataValue);
}

void precicef_read_bvdata_
(
  const int* dataID,
  const int* size,
  int*       valueIndices,
  double*    values )
{
  assertion(impl != NULL);
  impl->readBlockVectorData(*dataID, *size, valueIndices, values);
}

void precicef_read_vdata_
(
  const int* dataID,
  const int* valueIndex,
  double*    dataValue )
{
  assertion(impl != NULL);
  impl->readVectorData(*dataID, *valueIndex, dataValue);
}

void precicef_read_bsdata_
(
  const int* dataID,
  const int* size,
  int*       valueIndices,
  double*    values )
{
  assertion(impl != NULL);
  impl->readBlockScalarData(*dataID, *size, valueIndices, values);
}

void precicef_read_sdata_
(
  const int* dataID,
  const int* valueIndex,
  double*    dataValue )
{
  assertion(impl != NULL);
  impl->readScalarData(*dataID, *valueIndex, *dataValue);
}

void precicef_map_written_data_
(
  const int* meshID )
{
  assertion(impl != NULL);
  impl->mapWrittenData(*meshID);
}

void precicef_map_read_data_
(
  const int* meshID )
{
  assertion(impl != NULL);
  impl->mapReadData(*meshID);
}

void precicef_export_mesh_
(
  const char* filenameSuffix,
  int         filenameSuffixLength )
{
  assertion(impl != NULL);
  int strippedLength =
      precice::impl::strippedLength(filenameSuffix, filenameSuffixLength);
  string stringFilenameSuffix(filenameSuffix, strippedLength);
  impl->exportMesh(stringFilenameSuffix);
}

int precice::impl::strippedLength
(
  const char* string,
  int         length )
{
  int i=length-1;
  while (((string[i] == ' ') || (string[i] == 0)) && (i >= 0)){
    i--;
  }
  return i+1;
}

