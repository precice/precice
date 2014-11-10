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
  com::PtrCommunication com )
:
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
  const std::string& nameRequester,
  int                acceptorProcessRank,
  int                acceptorCommunicatorSize )
{
  preciceTrace2("acceptConnection()", nameAcceptor, nameRequester);
  if(not utils::MasterSlave::_slaveMode){
    if(utils::MasterSlave::_masterMode){
      _com->acceptConnection(nameAcceptor, nameRequester, acceptorProcessRank, 1);
    }
    else {
      _com->acceptConnection(nameAcceptor, nameRequester, acceptorProcessRank, acceptorCommunicatorSize);
    }
    _isConnected = _com->isConnected();
  }

  //broadcast isConnected
  if(utils::MasterSlave::_slaveMode){
    utils::MasterSlave::_communication->receive(_isConnected,0);
  }
  if(utils::MasterSlave::_masterMode){
    for(int rankSlave = 1; rankSlave < utils::MasterSlave::_size; rankSlave++){
      utils::MasterSlave::_communication->send(_isConnected,rankSlave);
    }
  }
}

void GatherScatterCommunication:: requestConnection
(
  const std::string& nameAcceptor,
  const std::string& nameRequester,
  int                requesterProcessRank,
  int                requesterCommunicatorSize )
{
  preciceTrace2("requestConnection()", nameAcceptor, nameRequester);
  if(not utils::MasterSlave::_slaveMode){
    if(utils::MasterSlave::_masterMode){
      _com->requestConnection(nameAcceptor, nameRequester, requesterProcessRank, 1);
    }
    else {
      _com->requestConnection(nameAcceptor, nameRequester, requesterProcessRank, requesterCommunicatorSize);
    }
    _isConnected = _com->isConnected();
  }

  //broadcast isConnected
  if(utils::MasterSlave::_slaveMode){
    utils::MasterSlave::_communication->receive(_isConnected,0);
  }
  if(utils::MasterSlave::_masterMode){
    for(int rankSlave = 1; rankSlave < utils::MasterSlave::_size; rankSlave++){
      utils::MasterSlave::_communication->send(_isConnected,rankSlave);
    }
  }
}

void GatherScatterCommunication:: closeConnection()
{
  preciceTrace("closeConnection()");
  if(not utils::MasterSlave::_slaveMode){
    _com->closeConnection();
  }
  //broadcast isConnected
  if(utils::MasterSlave::_slaveMode){
    utils::MasterSlave::_communication->receive(_isConnected,0);
  }
  if(utils::MasterSlave::_masterMode){
    for(int rankSlave = 1; rankSlave < utils::MasterSlave::_size; rankSlave++){
      utils::MasterSlave::_communication->send(_isConnected,rankSlave);
    }
  }
}

com::PtrCommunication GatherScatterCommunication:: getMasterCommunication()
{
  return _com;
}

void GatherScatterCommunication:: startSendPackage
(
  int rankReceiver )
{
  if(not utils::MasterSlave::_slaveMode){
    _com->startSendPackage(rankReceiver);
  }
}

void GatherScatterCommunication:: finishSendPackage()
{
  if(not utils::MasterSlave::_slaveMode){
   _com->finishSendPackage();
  }
}

int GatherScatterCommunication:: startReceivePackage
(
  int rankSender )
{
  preciceTrace1("startReceivePackage()", rankSender);
  if(not utils::MasterSlave::_slaveMode){
    return _com->startReceivePackage(rankSender);
  }
  return -1;
}

void GatherScatterCommunication:: finishReceivePackage()
{
  _com->finishReceivePackage();
}

void GatherScatterCommunication:: sendMaster
(
  const std::string& itemToSend,
  int                rankReceiver )
{
  preciceTrace2("send(string)", itemToSend, rankReceiver);
  if(not utils::MasterSlave::_slaveMode){
    _com->send(itemToSend, rankReceiver);
  }
}

void GatherScatterCommunication:: sendMaster
(
  int* itemsToSend,
  int  size,
  int  rankReceiver )
{
  preciceTrace2("send(int*)", size, rankReceiver);
  if(not utils::MasterSlave::_slaveMode){
    _com->send(itemsToSend, size, rankReceiver);
  }
}

void GatherScatterCommunication:: sendMaster
(
  double* itemsToSend,
  int     size,
  int     rankReceiver )
{
  preciceTrace2("send(double*)", size, rankReceiver);
  if(not utils::MasterSlave::_slaveMode){
    _com->send(itemsToSend, size, rankReceiver);
  }
}

void GatherScatterCommunication:: sendMaster
(
  double itemToSend,
  int    rankReceiver )
{
  preciceTrace2("send(double)", itemToSend, rankReceiver);
  if(not utils::MasterSlave::_slaveMode){
    _com->send(itemToSend, rankReceiver);
  }
}

void GatherScatterCommunication:: sendMaster
(
  int itemToSend,
  int rankReceiver )
{
  preciceTrace2("send(int)", itemToSend, rankReceiver);
  if(not utils::MasterSlave::_slaveMode){
    _com->send(itemToSend, rankReceiver);
  }
}

void GatherScatterCommunication:: sendMaster
(
  bool itemToSend,
  int  rankReceiver )
{
  preciceTrace2("send(bool)", itemToSend, rankReceiver);
  if(not utils::MasterSlave::_slaveMode){
    _com->send(itemToSend, rankReceiver);
  }
}

int GatherScatterCommunication:: receiveMaster
(
  std::string& itemToReceive,
  int          rankSender )
{
  preciceTrace1("receive(string)", rankSender);
  if(not utils::MasterSlave::_slaveMode){
    return _com->receive(itemToReceive, rankSender);
  }
  return -1;
}

int GatherScatterCommunication:: receiveMaster
(
  int* itemsToReceive,
  int  size,
  int  rankSender )
{
  preciceTrace2("receive(int*)", size, rankSender);
  if(not utils::MasterSlave::_slaveMode){
    return _com->receive(itemsToReceive, size, rankSender);
  }
  return -1;
}

int GatherScatterCommunication:: receiveMaster
(
  double* itemsToReceive,
  int     size,
  int     rankSender )
{
  preciceTrace2("receive(double*)", size, rankSender);
  if(not utils::MasterSlave::_slaveMode){
    return _com->receive(itemsToReceive, size, rankSender);
  }
  return -1;
}

int GatherScatterCommunication:: receiveMaster
(
  double& itemToReceive,
  int     rankSender )
{
  preciceTrace1("receive(double)", rankSender);
  if(not utils::MasterSlave::_slaveMode){
    return _com->receive(itemToReceive, rankSender);
  }
  return -1;
}

int GatherScatterCommunication:: receiveMaster
(
  int& itemToReceive,
  int  rankSender )
{
  preciceTrace1("receive(int)", rankSender);
  if(not utils::MasterSlave::_slaveMode){
    return _com->receive(itemToReceive, rankSender);
  }
  return -1;
}

int GatherScatterCommunication:: receiveMaster
(
  bool& itemToReceive,
  int   rankSender )
{
  preciceTrace1("receive(bool)", rankSender);
  if(not utils::MasterSlave::_slaveMode){
    return _com->receive(itemToReceive, rankSender);
  }
  return -1;
}


void GatherScatterCommunication:: sendAll (
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

      //master data
      for(int i=0; i<vertexDistribution[0].size();i++){
        for(int j=0;j<valueDimension;j++){
          globalItemsToSend[vertexDistribution[0][i]*valueDimension+j] += (*itemsToSend)[i*valueDimension+j];
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
void GatherScatterCommunication:: receiveAll (
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

      //master data
      for(int i=0; i<vertexDistribution[0].size();i++){
        for(int j=0;j<valueDimension;j++){
          (*itemsToReceive)[i*valueDimension+j] = globalItemsToReceive[vertexDistribution[0][i]*valueDimension+j];
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
          delete valuesSlave;
        }
      }
      delete globalItemsToReceive;
    } //master
  }
}

void GatherScatterCommunication:: receiveAll (
  bool&  itemToReceive,
  int    rankSender )
{
  preciceTrace1("receiveAll(bool)", utils::MasterSlave::_rank);
  if(not utils::MasterSlave::_slaveMode){
    _com->receive(itemToReceive, rankSender);
  }

  if(utils::MasterSlave::_slaveMode){
    assertion(utils::MasterSlave::_communication->isConnected());
    utils::MasterSlave::_communication->receive(itemToReceive,0);
  }
  if(utils::MasterSlave::_masterMode){
    assertion(utils::MasterSlave::_communication->isConnected());
    for(int rankSlave = 1; rankSlave < utils::MasterSlave::_size; rankSlave++){
      utils::MasterSlave::_communication->send(itemToReceive,rankSlave);
    }
  }
  preciceDebug("ReceiveAll(bool): " << itemToReceive);
}

void GatherScatterCommunication:: receiveAll (
  double&  itemToReceive,
  int      rankSender )
{
  preciceTrace1("receiveAll(double)", utils::MasterSlave::_rank);
  if(not utils::MasterSlave::_slaveMode){
    _com->receive(itemToReceive, rankSender);
  }

  if(utils::MasterSlave::_slaveMode){
    assertion(utils::MasterSlave::_communication->isConnected());
    utils::MasterSlave::_communication->receive(itemToReceive,0);
  }
  if(utils::MasterSlave::_masterMode){
    assertion(utils::MasterSlave::_communication->isConnected());
    for(int rankSlave = 1; rankSlave < utils::MasterSlave::_size; rankSlave++){
      utils::MasterSlave::_communication->send(itemToReceive,rankSlave);
    }
  }
  preciceDebug("ReceiveAll(bool): " << itemToReceive);
}

void GatherScatterCommunication:: sendAll (
  bool   itemToSend,
  int    rankReceiver)
{
  preciceTrace1("sendAll(bool)", utils::MasterSlave::_rank);
  if(not utils::MasterSlave::_slaveMode){
    _com->send(itemToSend, rankReceiver);
  }
}

void GatherScatterCommunication:: sendAll (
  double   itemToSend,
  int      rankReceiver)
{
  preciceTrace1("sendAll(double)", utils::MasterSlave::_rank);
  if(not utils::MasterSlave::_slaveMode){
    _com->send(itemToSend, rankReceiver);
  }
}

}} // namespace precice, m2n

