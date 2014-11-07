// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License

#include "SimpleCommunication.hpp"
#include "com/Communication.hpp"
#include "utils/MasterSlave.hpp"
#include "mesh/Mesh.hpp"

namespace precice {
namespace m2n {

tarch::logging::Log SimpleCommunication:: _log("precice::m2n::SimpleCommunication");

SimpleCommunication:: SimpleCommunication
(
  com::PtrCommunication com )
:
  _com(com)
{}

SimpleCommunication:: ~SimpleCommunication()
{
  if (isConnected()){
    closeConnection();
  }
}

bool SimpleCommunication:: isConnected()
{
  return _com->isConnected();
}

void SimpleCommunication:: acceptConnection
(
  const std::string& nameAcceptor,
  const std::string& nameRequester,
  int                acceptorProcessRank,
  int                acceptorCommunicatorSize )
{
  preciceTrace2("acceptConnection()", nameAcceptor, nameRequester);
  _com->acceptConnection(nameAcceptor, nameRequester, acceptorProcessRank, acceptorCommunicatorSize);
}

void SimpleCommunication:: requestConnection
(
  const std::string& nameAcceptor,
  const std::string& nameRequester,
  int                requesterProcessRank,
  int                requesterCommunicatorSize )
{
  preciceTrace2("requestConnection()", nameAcceptor, nameRequester);
  _com->requestConnection(nameAcceptor, nameRequester, requesterProcessRank, requesterCommunicatorSize);
}

void SimpleCommunication:: closeConnection()
{
  preciceTrace("closeConnection()");
  _com->closeConnection();
}

com::PtrCommunication SimpleCommunication:: getMasterCommunication()
{
  return _com;
}

void SimpleCommunication:: startSendPackage
(
  int rankReceiver )
{
  _com->startSendPackage(rankReceiver);
}

void SimpleCommunication:: finishSendPackage()
{
  _com->finishSendPackage();
}

int SimpleCommunication:: startReceivePackage
(
  int rankSender )
{
  preciceTrace1("startReceivePackage()", rankSender);
  return _com->startReceivePackage(rankSender);
}

void SimpleCommunication:: finishReceivePackage()
{
  _com->finishReceivePackage();
}

void SimpleCommunication:: sendMaster
(
  const std::string& itemToSend,
  int                rankReceiver )
{
  preciceTrace2("send(string)", itemToSend, rankReceiver);
  _com->send(itemToSend, rankReceiver);
}

void SimpleCommunication:: sendMaster
(
  int* itemsToSend,
  int  size,
  int  rankReceiver )
{
  preciceTrace2("send(int*)", size, rankReceiver);
  _com->send(itemsToSend, size, rankReceiver);
}

void SimpleCommunication:: sendMaster
(
  double* itemsToSend,
  int     size,
  int     rankReceiver )
{
  preciceTrace2("send(double*)", size, rankReceiver);
  _com->send(itemsToSend, size, rankReceiver);
}

void SimpleCommunication:: sendMaster
(
  double itemToSend,
  int    rankReceiver )
{
  preciceTrace2("send(double)", itemToSend, rankReceiver);
  _com->send(itemToSend, rankReceiver);
}

void SimpleCommunication:: sendMaster
(
  int itemToSend,
  int rankReceiver )
{
  preciceTrace2("send(int)", itemToSend, rankReceiver);
  _com->send(itemToSend, rankReceiver);
}

void SimpleCommunication:: sendMaster
(
  bool itemToSend,
  int  rankReceiver )
{
  preciceTrace2("send(bool)", itemToSend, rankReceiver);
  _com->send(itemToSend, rankReceiver);
}

int SimpleCommunication:: receiveMaster
(
  std::string& itemToReceive,
  int          rankSender )
{
  preciceTrace1("receive(string)", rankSender);
  return _com->receive(itemToReceive, rankSender);
}

int SimpleCommunication:: receiveMaster
(
  int* itemsToReceive,
  int  size,
  int  rankSender )
{
  preciceTrace2("receive(int*)", size, rankSender);
  return _com->receive(itemsToReceive, size, rankSender);
}

int SimpleCommunication:: receiveMaster
(
  double* itemsToReceive,
  int     size,
  int     rankSender )
{
  preciceTrace2("receive(double*)", size, rankSender);
  return _com->receive(itemsToReceive, size, rankSender);
}

int SimpleCommunication:: receiveMaster
(
  double& itemToReceive,
  int     rankSender )
{
  preciceTrace1("receive(double)", rankSender);
  return _com->receive(itemToReceive, rankSender);
}

int SimpleCommunication:: receiveMaster
(
  int& itemToReceive,
  int  rankSender )
{
  preciceTrace1("receive(int)", rankSender);
  return _com->receive(itemToReceive, rankSender);
}

int SimpleCommunication:: receiveMaster
(
  bool& itemToReceive,
  int   rankSender )
{
  preciceTrace1("receive(bool)", rankSender);
  return _com->receive(itemToReceive, rankSender);
}


void SimpleCommunication:: sendAll (
  utils::DynVector*    itemsToSend,
  int           size,
  int           rankReceiver,
  mesh::PtrMesh mesh,
  int           valueDimension)
{
  preciceTrace1("sendAll", size);

  double* globalItemsToSend = NULL;
  int globalSize = -1;

  //gatherData
  if(utils::MasterSlave::_slaveMode || utils::MasterSlave::_masterMode){

    assertion(utils::MasterSlave::_communication.get() != NULL);
    assertion(utils::MasterSlave::_communication->isConnected());
    assertion(utils::MasterSlave::_size>1);
    assertion(utils::MasterSlave::_rank!=-1);

    if(utils::MasterSlave::_rank>0){ //slave
      if (size > 0) {
        utils::MasterSlave::_communication->send(tarch::la::raw(*itemsToSend), size, 0);
      }
    }
    else{ //master
      assertion(utils::MasterSlave::_rank==0);
      std::map<int,std::vector<int> >& vertexDistribution = mesh->getVertexDistribution();
      globalSize = mesh->getGlobalNumberOfVertices()*valueDimension;
      preciceDebug("Global Size = " << globalSize);
      globalItemsToSend = new double[globalSize]();

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
          delete valuesSlave;
        }
      }
    } //master
  }
  else{ //couplingMode
    globalItemsToSend = tarch::la::raw(*itemsToSend);
    globalSize = size;
  }

  //send Data to other participant
  if(not utils::MasterSlave::_slaveMode ){
    assertion(globalItemsToSend!=NULL);
    _com->send(globalItemsToSend, globalSize, rankReceiver);
  }

  if(utils::MasterSlave::_masterMode){
    delete globalItemsToSend;
  }
}


/**
 * @brief Receives an array of double values.
 *
 * @return Rank of sender, which is useful when ANY_SENDER is used.
 */
void SimpleCommunication:: receiveAll (
  utils::DynVector*   itemsToReceive,
  int           size,
  int           rankSender,
  mesh::PtrMesh mesh,
  int           valueDimension)
{
  preciceTrace1("receiveAll", size);

  double* globalItemsToReceive = NULL;
  int globalSize = -1;

  // receive Data from other participant
  if(not utils::MasterSlave::_slaveMode ){
    if(utils::MasterSlave::_masterMode){
      globalSize = mesh->getGlobalNumberOfVertices()*valueDimension;
      preciceDebug("Global Size = " << globalSize);
      globalItemsToReceive = new double[globalSize];
    }
    else{ //couplingMode
      globalItemsToReceive = tarch::la::raw(*itemsToReceive);
      globalSize = size;
    }
    assertion(globalItemsToReceive!=NULL);
    _com->receive(globalItemsToReceive, globalSize, rankSender);
  }

  // scatter Data
  if(utils::MasterSlave::_slaveMode || utils::MasterSlave::_masterMode){

    assertion(utils::MasterSlave::_communication.get() != NULL);
    assertion(utils::MasterSlave::_communication->isConnected());
    assertion(utils::MasterSlave::_size>1);
    assertion(utils::MasterSlave::_rank!=-1);

    if(utils::MasterSlave::_rank>0){ //slave
      if (size > 0) {
        utils::MasterSlave::_communication->receive(tarch::la::raw(*itemsToReceive), size, 0);
      }
    }
    else{ //master
      assertion(utils::MasterSlave::_rank==0);
      std::map<int,std::vector<int> >& vertexDistribution = mesh->getVertexDistribution();

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
          delete valuesSlave;
        }
      }
      delete globalItemsToReceive;
    } //master
  }
}

}} // namespace precice, m2n

