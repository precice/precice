// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License

#include "GatherScatterCommunication.hpp"
#include "com/Communication.hpp"
#include "utils/MasterSlave.hpp"
#include "mesh/Mesh.hpp"

namespace precice {
namespace m2n {

tarch::logging::Log GatherScatterCommunication:: _log("precice::m2n::GatherScatterCommunication");

GatherScatterCommunication:: GatherScatterCommunication
(
  com::Communication::SharedPointer com,
  mesh::PtrMesh mesh)
:
  DistributedCommunication(mesh),
  _com(com),
  _isConnected(false)
{}

GatherScatterCommunication:: ~GatherScatterCommunication()
{
  if (isConnected()){
    closeConnection();
  }
}

bool GatherScatterCommunication:: isConnected()
{
  return _isConnected;
}

void GatherScatterCommunication:: acceptConnection
(
  const std::string& nameAcceptor,
  const std::string& nameRequester)
{
  preciceTrace2("acceptConnection()", nameAcceptor, nameRequester);
  assertion(utils::MasterSlave::_slaveMode || _com->isConnected());
  _isConnected = true;
}

void GatherScatterCommunication:: requestConnection
(
  const std::string& nameAcceptor,
  const std::string& nameRequester )
{
  preciceTrace2("requestConnection()", nameAcceptor, nameRequester);
  assertion(utils::MasterSlave::_slaveMode || _com->isConnected());
  _isConnected = true;
}

void GatherScatterCommunication:: closeConnection()
{
  preciceTrace("closeConnection()");
  assertion(utils::MasterSlave::_slaveMode || not _com->isConnected());
  _isConnected = false;
}


void GatherScatterCommunication:: send (
  double*    itemsToSend,
  int        size,
  int        valueDimension)
{
  preciceTrace1("sendAll", size);
  assertion(utils::MasterSlave::_slaveMode || utils::MasterSlave::_masterMode);
  assertion(utils::MasterSlave::_communication.get() != NULL);
  assertion(utils::MasterSlave::_communication->isConnected());
  assertion(utils::MasterSlave::_size>1);
  assertion(utils::MasterSlave::_rank!=-1);

  double* globalItemsToSend = NULL;

  //gatherData
  if(utils::MasterSlave::_slaveMode){ //slave
    if (size > 0) {
      utils::MasterSlave::_communication->send(itemsToSend, size, 0);
    }
  }
  else{ //master
    assertion(utils::MasterSlave::_rank==0);
    std::map<int,std::vector<int> >& vertexDistribution = _mesh->getVertexDistribution();
    int globalSize = _mesh->getGlobalNumberOfVertices()*valueDimension;
    preciceDebug("Global Size = " << globalSize);
    globalItemsToSend = new double[globalSize]();

    //master data
    for(int i=0; i<vertexDistribution[0].size();i++){
      for(int j=0;j<valueDimension;j++){
        globalItemsToSend[vertexDistribution[0][i]*valueDimension+j] += itemsToSend[i*valueDimension+j];
      }
    }

    //slaves data
    for(int rankSlave = 1; rankSlave < utils::MasterSlave::_size; rankSlave++){
      int slaveSize = vertexDistribution[rankSlave].size()*valueDimension;
      preciceDebug("Slave Size = " << slaveSize );
      if (slaveSize > 0) {
        double* valuesSlave = new double[slaveSize];
        utils::MasterSlave::_communication->receive(valuesSlave, slaveSize, rankSlave);
        for(int i=0; i<vertexDistribution[rankSlave].size();i++){
          for(int j=0;j<valueDimension;j++){
            globalItemsToSend[vertexDistribution[rankSlave][i]*valueDimension+j] += valuesSlave[i*valueDimension+j];
          }
        }
        delete[] valuesSlave;
      }
    }

    //send data to other master
    assertion(globalItemsToSend!=NULL);
    _com->send(globalItemsToSend, globalSize, 0);
    delete[] globalItemsToSend;
  } //master
}


/**
 * @brief Receives an array of double values.
 *
 * @return Rank of sender, which is useful when ANY_SENDER is used.
 */
void GatherScatterCommunication:: receive (
  double*   itemsToReceive,
  int       size,
  int       valueDimension)
{
  preciceTrace1("receiveAll", size);
  assertion(utils::MasterSlave::_slaveMode || utils::MasterSlave::_masterMode);
  assertion(utils::MasterSlave::_communication.get() != NULL);
  assertion(utils::MasterSlave::_communication->isConnected());
  assertion(utils::MasterSlave::_size>1);
  assertion(utils::MasterSlave::_rank!=-1);

  double* globalItemsToReceive = NULL;

  //receive data at master
  if(utils::MasterSlave::_masterMode){
    int globalSize = _mesh->getGlobalNumberOfVertices()*valueDimension;
    preciceDebug("Global Size = " << globalSize);
    globalItemsToReceive = new double[globalSize];
    _com->receive(globalItemsToReceive, globalSize, 0);
  }

  //scatter data
  if(utils::MasterSlave::_slaveMode){ //slave
    if (size > 0) {
      preciceDebug("itemsToRec[0] = " << itemsToReceive[0]);
      utils::MasterSlave::_communication->receive(itemsToReceive, size, 0);
      preciceDebug("itemsToRec[0] = " << itemsToReceive[0]);
    }
  }
  else{ //master
    assertion(utils::MasterSlave::_rank==0);
    std::map<int,std::vector<int> >& vertexDistribution = _mesh->getVertexDistribution();

    //master data
    for(int i=0; i<vertexDistribution[0].size();i++){
      for(int j=0;j<valueDimension;j++){
        itemsToReceive[i*valueDimension+j] = globalItemsToReceive[vertexDistribution[0][i]*valueDimension+j];
      }
    }

    //slaves data
    for(int rankSlave = 1; rankSlave < utils::MasterSlave::_size; rankSlave++){
      int slaveSize = vertexDistribution[rankSlave].size()*valueDimension;
      preciceDebug("Slave Size = " << slaveSize );
      if (slaveSize > 0) {
        double* valuesSlave = new double[slaveSize];
        for(int i=0; i<vertexDistribution[rankSlave].size();i++){
          for(int j=0;j<valueDimension;j++){
            valuesSlave[i*valueDimension+j] = globalItemsToReceive[vertexDistribution[rankSlave][i]*valueDimension+j];
          }
        }
        utils::MasterSlave::_communication->send(valuesSlave, slaveSize, rankSlave);
        preciceDebug("valuesSlave[0] = " << valuesSlave[0]);
        delete[] valuesSlave;
      }
    }
    delete[] globalItemsToReceive;
  } //master
}

}} // namespace precice, m2n
