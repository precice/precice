// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_NO_MPI

#include "MPIPortsCommunication.hpp"
#include "utils/Globals.hpp"
#include "utils/Parallel.hpp"
#include <fstream>

namespace precice {
namespace com {

tarch::logging::Log MPIPortsCommunication:: _log ("precice::com::MPIPortsCommunication");


MPIPortsCommunication:: MPIPortsCommunication
(
  const std::string& publishingDirectory )
:
  MPICommunication ( utils::Parallel::getGlobalCommunicator() ),
  _publishingDirectory(publishingDirectory),
  _portname (),
  _isConnection ( false )
{}

MPIPortsCommunication:: ~MPIPortsCommunication()
{
  preciceTrace1 ( "~MPIPortsCommunication()", _isConnection );
  if (_isConnection){
    closeConnection();
  }
}

int MPIPortsCommunication:: getRemoteCommunicatorSize()
{
  preciceTrace ( "getRemoteCommunicatorSize()" );
  assertion ( _isConnection );
  int remoteSize = 0;
  MPI_Comm_remote_size ( communicator(), &remoteSize );
  return remoteSize;
}

void MPIPortsCommunication:: acceptConnection
(
  const std::string& nameAcceptor,
  const std::string& nameRequester,
  int                acceptorProcessRank,
  int                acceptorCommunicatorSize )
{
  preciceTrace2 ( "acceptConnection()", nameAcceptor, nameRequester );
  assertion ( not _isConnection );

  int argc = 1;
  char* arg = new char[8];
  strcpy(arg, "precice");
  char** argv = &arg;
  utils::Parallel::initialize(&argc, &argv, nameAcceptor);
  delete[] arg;

  MPI_Open_port(MPI_INFO_NULL, _portname);
  // Write portname to file
  std::string portFilename ( _publishingDirectory + "." + nameRequester + "-portname" );
  preciceDebug ( "Writing server connection info to file " + portFilename );
  std::ofstream outFile;
  outFile.open ((portFilename +  "~").c_str(), std::ios::out);
  outFile << _portname;
  outFile.close();
  // To give the file first a "wrong" name prevents early reading errors
  rename( (portFilename + "~").c_str(), portFilename.c_str() );

  preciceDebug("Calling MPI_Comm_accept() with portname = " << _portname);
  MPI_Comm localComm = utils::Parallel::getLocalCommunicator();
  MPI_Comm_accept ( _portname, MPI_INFO_NULL, 0, localComm, &communicator() );
  if ( utils::Parallel::getLocalProcessRank() == 0 ){
    if ( remove(portFilename.c_str()) != 0 ) {
      preciceWarning ( "acceptConnection()", "Could not remove port information file!" );
    }
  }
  _isConnection = true;
}

void MPIPortsCommunication:: requestConnection
(
  const std::string& nameAcceptor,
  const std::string& nameRequester,
  int                requesterProcessRank,
  int                requesterCommunicatorSize )
{
  preciceTrace2 ( "requestConnection()", nameAcceptor, nameRequester );
  assertion ( not _isConnection );

  int argc = 1;
  char* arg = new char[8];
  strcpy(arg, "precice");
  char** argv = &arg;
  utils::Parallel::initialize ( &argc, &argv, nameRequester );
  delete[] arg;

  std::string portFilename(_publishingDirectory + "." + nameRequester + "-portname");
  std::ifstream inFile;
  do {
    inFile.open ( portFilename.c_str(), std::ios::in );
  } while ( not inFile );
  inFile.getline ( _portname, MPI_MAX_PORT_NAME );
  inFile.close();
  preciceDebug ( "Read connection info \"" + std::string(_portname) +
                 "\" from file " + portFilename );
  preciceDebug("Calling MPI_Comm_connect() with portname = " << _portname);
  MPI_Comm localComm = utils::Parallel::getLocalCommunicator();
  MPI_Comm_connect ( _portname, MPI_INFO_NULL, 0, localComm, &communicator() );
  _isConnection = true;
}

void MPIPortsCommunication:: closeConnection()
{
  assertion ( _isConnection );
  MPI_Comm_disconnect ( & communicator() );
  _isConnection = false;
}

}} // namespace precice, com

#endif // not PRECICE_NO_MPI
