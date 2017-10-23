#include "SolverInterfaceFASTEST.hpp"
#include "SolverInterfaceFortran.hpp"
#include "precice/impl/SolverInterfaceImpl.hpp"
#include <iostream>
#include <string>
#include "logging/Logger.hpp"

using namespace std;

static precice::impl::SolverInterfaceImpl* implA = nullptr;
static precice::impl::SolverInterfaceImpl* implF = nullptr;

static precice::logging::Logger _log ("SolverInterfaceFASTEST");

static std::string errormsg = "preCICE has not been created properly. Be sure to call \"precicef_create\" before any other call to preCICE.";

namespace precice {
  namespace impl {
    /**
     * @brief Returns length of string without trailing whitespaces.
     */
    int strippedLength( const char* string, int length );
  }
}

void precice_fastest_create_
(
  const char* accessorNameA,
  const char* accessorNameF,
  const char* configFileName,
  const int*  solverProcessIndex,
  const int*  solverProcessSize,
  int   lengthAccessorNameA,
  int   lengthAccessorNameF,
  int   lengthConfigFileName )
{
  int strippedLength = precice::impl::strippedLength(accessorNameA,lengthAccessorNameA);
  string stringAccessorNameA(accessorNameA, strippedLength);
  strippedLength = precice::impl::strippedLength(accessorNameF,lengthAccessorNameF);
  string stringAccessorNameF(accessorNameF, strippedLength);
  strippedLength = precice::impl::strippedLength(configFileName,lengthConfigFileName);
  string stringConfigFileName(configFileName, strippedLength);
  //cout << "Accessor: " << stringAccessorName << "!" << endl;
  //cout << "Config  : " << stringConfigFileName << "!" << endl;
  implA = new precice::impl::SolverInterfaceImpl (stringAccessorNameA,
         *solverProcessIndex, *solverProcessSize, false);
  implA->configure(stringConfigFileName);
  implF = new precice::impl::SolverInterfaceImpl (stringAccessorNameF,
         *solverProcessIndex, *solverProcessSize, false);
  implF->configure(stringConfigFileName);
}

void precice_fastest_initialize_
(
  double* timestepLengthLimit,
  const int*  useF)
{
  CHECK(implA != nullptr && implF != nullptr ,errormsg);
  assertion(*useF == 0 || *useF == 1);
  if(*useF==0){
    *timestepLengthLimit = implA->initialize();
  }
  else{
    *timestepLengthLimit = implF->initialize();
  }
}

void precice_fastest_initialize_data_(
  const int*  useF)
{
  CHECK(implA != nullptr && implF != nullptr ,errormsg);
  assertion(*useF == 0 || *useF == 1);
  if(*useF==0){
    implA->initializeData();
  }
  else{
    implF->initializeData();
  }
}

void precice_fastest_advance_
(
  double* timestepLengthLimit,
  const int*  useF)
{
  CHECK(implA != nullptr && implF != nullptr ,errormsg);
  assertion(*useF == 0 || *useF == 1);
  if(*useF==0){
    *timestepLengthLimit = implA->advance(*timestepLengthLimit);
  }
  else{
    *timestepLengthLimit = implF->advance(*timestepLengthLimit);
  }
}

void precice_fastest_finalize_(
  const int*  useF)
{
  CHECK(implA != nullptr && implF != nullptr ,errormsg);
  assertion(*useF == 0 || *useF == 1);
  if(*useF==0){
    implA->finalize();
    delete implA;
  }
  else{
    implF->finalize();
    delete implF;
  }
}

void precice_fastest_action_required_
(
  const char* action,
  int*        isRequired,
  const int*  useF,
  int         lengthAction )
{
  CHECK(implA != nullptr && implF != nullptr ,errormsg);
  assertion(*useF == 0 || *useF == 1);

  int strippedLength = precice::impl::strippedLength(action, lengthAction);
  string stringAction(action, strippedLength);
  if(*useF==0){
    if (implA->isActionRequired(stringAction)){
      *isRequired = 1;
    }
    else {
      *isRequired = 0;
    }
  }
  else{
    if (implF->isActionRequired(stringAction)){
      *isRequired = 1;
    }
    else {
      *isRequired = 0;
    }
  }
}

void precice_fastest_fulfilled_action_
(
  const char* action,
  const int*  useF,
  int         lengthAction )
{
  CHECK(implA != nullptr && implF != nullptr ,errormsg);
  assertion(*useF == 0 || *useF == 1);
  int strippedLength = precice::impl::strippedLength(action, lengthAction);
  string stringAction(action, strippedLength);
  if(*useF==0){
    implA->fulfilledAction(stringAction);
  }
  else{
    implF->fulfilledAction(stringAction);
  }
}

void precice_fastest_get_mesh_id_
(
  const char* meshName,
  int*        meshID,
  const int*  useF,
  int         lengthMeshName )
{
  CHECK(implA != nullptr && implF != nullptr ,errormsg);
  assertion(*useF == 0 || *useF == 1);
  int strippedLength = precice::impl::strippedLength(meshName, lengthMeshName);
  string stringMeshName(meshName, strippedLength);
  if(*useF==0){
    *meshID = implA->getMeshID(stringMeshName);
  }
  else{
    *meshID = implF->getMeshID(stringMeshName);
  }
}

void precice_fastest_get_data_id_
(
  const char* dataName,
  const int*  meshID,
  int*        dataID,
  const int*  useF,
  int         lengthDataName
)
{
  CHECK(implA != nullptr && implF != nullptr ,errormsg);
  assertion(*useF == 0 || *useF == 1);
  int strippedLength = precice::impl::strippedLength(dataName, lengthDataName);
  string stringDataName(dataName, strippedLength);
  if(*useF==0){
    *dataID = implA->getDataID(stringDataName, *meshID);
  }
  else{
    *dataID = implF->getDataID(stringDataName, *meshID);
  }
}

void precice_fastest_set_vertex_
(
  const int*    meshID,
  const double* position,
  int*          vertexID,
  const int*    useF)
{
  CHECK(implA != nullptr && implF != nullptr ,errormsg);
  assertion(*useF == 0 || *useF == 1);
  if(*useF==0){
    *vertexID = implA->setMeshVertex(*meshID, position);
  }
  else{
    *vertexID = implF->setMeshVertex(*meshID, position);
  }
}

void precice_fastest_set_vertices_
(
  const int*    meshID,
  const int*    size,
  double*       positions,
  int*          positionIDs,
  const int*    useF)
{
  CHECK(implA != nullptr && implF != nullptr ,errormsg);
  assertion(*useF == 0 || *useF == 1);
  if(*useF==0){
    implA->setMeshVertices(*meshID, *size, positions, positionIDs);
  }
  else{
    implF->setMeshVertices(*meshID, *size, positions, positionIDs);
  }
}

void precice_fastest_write_bvdata_
(
  const int* dataID,
  const int* size,
  int*       valueIndices,
  double*    values,
  const int* useF)
{
  CHECK(implA != nullptr && implF != nullptr ,errormsg);
  assertion(*useF == 0 || *useF == 1);
  if(*useF==0){
    implA->writeBlockVectorData(*dataID, *size, valueIndices, values);
  }
  else{
    implF->writeBlockVectorData(*dataID, *size, valueIndices, values);
  }
}

void precice_fastest_write_vdata_
(
  const int*    dataID,
  const int*    valueIndex,
  const double* dataValue,
  const int*    useF)
{
  CHECK(implA != nullptr && implF != nullptr ,errormsg);
  assertion(*useF == 0 || *useF == 1);
  if(*useF==0){
    implA->writeVectorData(*dataID, *valueIndex, dataValue);
  }
  else{
    implF->writeVectorData(*dataID, *valueIndex, dataValue);
  }
}

void precice_fastest_write_bsdata_
(
  const int* dataID,
  const int* size,
  int*       valueIndices,
  double*    values,
  const int* useF)
{
  CHECK(implA != nullptr && implF != nullptr ,errormsg);
  assertion(*useF == 0 || *useF == 1);
  if(*useF==0){
    implA->writeBlockScalarData(*dataID, *size, valueIndices, values);
  }
  else{
    implF->writeBlockScalarData(*dataID, *size, valueIndices, values);
  }
}

void precice_fastest_write_sdata_
(
  const int*    dataID,
  const int*    valueIndex,
  const double* dataValue,
  const int*    useF)
{
  CHECK(implA != nullptr && implF != nullptr ,errormsg);
  assertion(*useF == 0 || *useF == 1);
  if(*useF==0){
    implA->writeScalarData(*dataID, *valueIndex, *dataValue);
  }
  else{
    implF->writeScalarData(*dataID, *valueIndex, *dataValue);
  }
}

void precice_fastest_read_bvdata_
(
  const int* dataID,
  const int* size,
  int*       valueIndices,
  double*    values,
  const int* useF)
{
  CHECK(implA != nullptr && implF != nullptr ,errormsg);
  assertion(*useF == 0 || *useF == 1);
  if(*useF==0){
    implA->readBlockVectorData(*dataID, *size, valueIndices, values);
  }
  else{
    implF->readBlockVectorData(*dataID, *size, valueIndices, values);
  }
}

void precice_fastest_read_vdata_
(
  const int* dataID,
  const int* valueIndex,
  double*    dataValue,
  const int* useF)
{
  CHECK(implA != nullptr && implF != nullptr ,errormsg);
  assertion(*useF == 0 || *useF == 1);
  if(*useF==0){
    implA->readVectorData(*dataID, *valueIndex, dataValue);
  }
  else{
    implF->readVectorData(*dataID, *valueIndex, dataValue);
  }
}

void precice_fastest_read_bsdata_
(
  const int* dataID,
  const int* size,
  int*       valueIndices,
  double*    values,
  const int* useF)
{
  CHECK(implA != nullptr && implF != nullptr ,errormsg);
  assertion(*useF == 0 || *useF == 1);
  if(*useF==0){
    implA->readBlockScalarData(*dataID, *size, valueIndices, values);
  }
  else{
    implF->readBlockScalarData(*dataID, *size, valueIndices, values);
  }
}

void precice_fastest_read_sdata_
(
  const int* dataID,
  const int* valueIndex,
  double*    dataValue,
  const int* useF)
{
  CHECK(implA != nullptr && implF != nullptr ,errormsg);
  assertion(*useF == 0 || *useF == 1);
  if(*useF==0){
    implA->readScalarData(*dataID, *valueIndex, *dataValue);
  }
  else{
    implF->readScalarData(*dataID, *valueIndex, *dataValue);
  }
}

//int precice::impl::strippedLength
//(
//  const char* string,
//  int         length )
//{
//  int i=length-1;
//  while (((string[i] == ' ') || (string[i] == 0)) && (i >= 0)){
//    i--;
//  }
//  return i+1;
//}

