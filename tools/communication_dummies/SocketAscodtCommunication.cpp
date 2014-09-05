
#include "SocketAscodtCommunication.h"

SocketAscodtCommunication::SocketAscodtCommunication() {
}

SocketAscodtCommunication::~SocketAscodtCommunication() {
}



bool SocketAscodtCommunication::isConnected()
{
  return false;
}


void SocketAscodtCommunication::acceptConnection (
    int* vertexTable,
    int* adressTable)
{
  //TODO include here also some handshake test if the connection between the 2 masters works
}

void SocketAscodtCommunication::requestConnection (
    int* vertexTable,
    int* adressTable)
{
  //TODO include here also some handshake test if the connection between the 2 masters works
}

void SocketAscodtCommunication::closeConnection()
{

}

void SocketAscodtCommunication::send (
  double* itemsToSend,
  int     size)
{

}

void SocketAscodtCommunication::receive (
  double* itemsToReceive,
  int     size)
{

}



