// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
//#ifndef PRECICE_NO_MPI

#include "Parallel.hpp"
#include "utils/Globals.hpp"
#include "com/MPIDirectCommunication.hpp"
#include <map>
#ifndef PRECICE_NO_PETSC
#include "petsc.h"
#endif

namespace precice {
namespace utils {

tarch::logging::Log Parallel:: _log ( "precice::utils::Parallel" );

Parallel::Communicator Parallel:: _globalCommunicator = Parallel::getCommunicatorWorld();

Parallel::Communicator Parallel:: _localCommunicator = MPI_COMM_NULL;

//bool Parallel:: _isLocalCommunicatorSet = false;

bool Parallel:: _isInitialized = false;

bool Parallel:: _autoInitialized = false;

std::vector<Parallel::AccessorGroup> Parallel:: _accessorGroups;

Parallel::Communicator Parallel:: getCommunicatorWorld()
{
# ifndef PRECICE_NO_MPI
  return MPI_COMM_WORLD;
# else
  return -1;
# endif
}

void Parallel:: initialize
(
  int*               argc,
  char***            argv,
  const std::string& groupName )
{
# ifndef PRECICE_NO_MPI
  preciceTrace1 ( "initialize()", groupName );
  int isMPIInitialized;
  MPI_Initialized ( & isMPIInitialized );
  if ( not isMPIInitialized ){
    preciceDebug ( "Initialize MPI" );
    _autoInitialized = true;
    MPI_Init ( argc, argv );
  }

  // Exchange group information
  if ( _accessorGroups.empty() ){
    preciceDebug ( "Exchange group information" );
    _accessorGroups.clear(); // Makes reinitialization possible
    std::map<std::string,int> groupMap; // map from names to group ID
    MPI_Comm globalComm = getGlobalCommunicator();
    int rank = -1;
    MPI_Comm_rank ( globalComm, &rank );
    int size = -1;
    MPI_Comm_size ( globalComm, &size );

    bool severalGroups = false;

    if ( size > 1 ){
      com::MPIDirectCommunication com;
      if ( rank == 0 ){
        groupMap[groupName] = 0;
        AccessorGroup newGroup;
        newGroup.id = 0;
        newGroup.size = 1;
        newGroup.leaderRank = 0;
        newGroup.name = groupName;
        _accessorGroups.push_back(newGroup);
        for ( int i = 1; i < size; i++ ){
          std::string name;
          com.receive ( name, i );
          if ( groupMap.find(name) == groupMap.end() ){
            groupMap[name] = _accessorGroups.size();
            AccessorGroup newGroup;
            newGroup.id = _accessorGroups.size();
            newGroup.size = 1;
            newGroup.leaderRank = i;
            newGroup.name = name;
            _accessorGroups.push_back(newGroup);
          }
          else {
            _accessorGroups[groupMap[name]].size++;
          }
        }
        int groupCount = (int)_accessorGroups.size();
        MPI_Bcast ( &groupCount, 1, MPI_INT, 0, globalComm );

        typedef std::map<std::string,int>::value_type Pair;
        foreach ( AccessorGroup& group, _accessorGroups ){
          for ( int i = 1; i < size; i++ ){
            com.send ( group.name, i );
            com.send ( group.leaderRank, i );
            com.send ( group.id, i );
            com.send ( group.size, i );
          }
        }
        severalGroups = _accessorGroups.size() > 1;
      }
      else {
        com.send ( groupName, 0 );
        int groupCount = -1;
        MPI_Bcast ( &groupCount, 1, MPI_INT, 0, globalComm );
        severalGroups = groupCount > 1;
        for ( int i=0; i < groupCount; i++ ){
          AccessorGroup newGroup;
          com.receive ( newGroup.name, 0 );
          com.receive ( newGroup.leaderRank, 0 );
          com.receive ( newGroup.id, 0 );
          com.receive ( newGroup.size, 0 );
          _accessorGroups.push_back ( newGroup );
        }
      }

      if ( severalGroups ){
        preciceDebug ( "Split groups" );
        int groupID = -1;
         for ( size_t i=0; i < _accessorGroups.size(); i++ ){
           if ( _accessorGroups[i].name == groupName ){
             groupID = _accessorGroups[i].id;
           }
        }
        assertion ( groupID != -1 );
        MPI_Comm_split ( _globalCommunicator, groupID,
                         getProcessRank(), & _localCommunicator );
        //_isLocalCommunicatorSet = true;
      }

#     ifdef Debug
      preciceDebug ( "Detected " << _accessorGroups.size() << " groups" );
      foreach ( const AccessorGroup& group, _accessorGroups ) {
        preciceDebug ( "Group " << group.id << ": name = " << group.name
                       << ", leaderRank = " << group.leaderRank
                       << ", size = " << group.size );
      }
#     endif // Debug
    }
  }
# endif // not PRECICE_NO_MPI
# ifndef PRECICE_NO_PETSC
  PetscInitialize(argc, argv, "", NULL);
# endif
  _isInitialized = true;
}

void Parallel:: finalize()
{
  //assertion(_isInitialized);
  _accessorGroups.clear();
# ifndef PRECICE_NO_MPI
  if (_autoInitialized){
    MPI_Finalize();
  }
# endif // not PRECICE_NO_MPI
# ifndef PRECICE_NO_PETSC
  PetscFinalize();
# endif
}

int Parallel:: getProcessRank()
{
  // Do not use preciceTrace or preciceDebug here!
  int processRank = 0;
# ifndef PRECICE_NO_MPI
  if ( _isInitialized ){
    MPI_Comm_rank (_globalCommunicator, &processRank);
  }
# endif // not PRECICE_NO_MPI
  return processRank;
}

int Parallel:: getLocalProcessRank()
{
  int processRank = 0;
# ifndef PRECICE_NO_MPI
  if ( _accessorGroups.size() > 1 ){
    //assertion(_isLocalCommunicatorSet);
    MPI_Comm_rank ( _localCommunicator, &processRank );
  }
  else {
    processRank = getProcessRank();
  }
# endif
  return processRank;
}

int Parallel:: getCommunicatorSize()
{
  preciceTrace ( "getCommunicatorSize()" );
  int communicatorSize = 1;
# ifndef PRECICE_NO_MPI
  if ( _isInitialized ){
    MPI_Comm_size (_globalCommunicator, &communicatorSize);
  }
# endif // not PRECICE_NO_MPI
  return communicatorSize;
}

void Parallel:: synchronizeProcesses()
{
# ifndef PRECICE_NO_MPI
  preciceTrace ( "synchronizeProcesses()" );
  assertion ( _isInitialized );
  MPI_Barrier ( _globalCommunicator );
# endif // not PRECICE_NO_MPI
}

void Parallel:: synchronizeLocalProcesses()
{
# ifndef PRECICE_NO_MPI
  preciceTrace ( "synchronizeLocalProcesses()" );
  if ( _accessorGroups.size() > 1 ){
    //assertion(_isLocalCommunicatorSet);
    MPI_Barrier ( _localCommunicator );
  }
  else {
    synchronizeProcesses();
  }
# endif // not PRECICE_NO_MPI
}

void Parallel:: setGlobalCommunicator
(
  Parallel::Communicator defaultCommunicator )
{
# ifndef PRECICE_NO_MPI
  preciceTrace ( "setGlobalCommunicator()" );
  if (_globalCommunicator != getCommunicatorWorld()){
    MPI_Comm_free ( & _globalCommunicator );
  }
  _globalCommunicator = defaultCommunicator;
  _localCommunicator = MPI_COMM_NULL;
  //_isLocalCommunicatorSet = false;
  _accessorGroups.clear();
# endif // not PRECICE_NO_MPI
}

const Parallel::Communicator& Parallel:: getGlobalCommunicator()
{
   preciceTrace ( "getGlobalCommunicator()" );
   return _globalCommunicator;
}

const Parallel::Communicator& Parallel:: getLocalCommunicator()
{
  preciceTrace ( "getLocalCommunicator()" );
  if (_localCommunicator == MPI_COMM_NULL) {
    return _globalCommunicator; // global == local, if no local is set
  }
  return _localCommunicator;
}

Parallel::Communicator Parallel:: getRestrictedCommunicator
(
  const std::vector<int>& ranks )
{
  preciceTrace ( "getRestrictedCommunicator()" );
  Communicator restrictedCommunicator = getCommunicatorWorld();
# ifndef PRECICE_NO_MPI
  assertion ( _isInitialized );
  assertion ( ranks.size() > 0 );
  // Create group, containing all processes of communicator
  MPI_Group currentGroup;
  MPI_Comm_group ( _globalCommunicator, &currentGroup );
  int * ranksArray = new int[ranks.size()];
# ifdef Debug
  int communicatorSize = 0;
  MPI_Comm_size (_globalCommunicator, &communicatorSize);
  preciceDebug ( "Comm size:" << communicatorSize );
  preciceDebug ( "ranks size:" << (int)ranks.size() );
  //   assertion ( (int)ranks.size() < communicatorSize );
# endif // Debug
  for ( size_t i=0; i < ranks.size(); i++ ) {
    preciceDebug ( "Adding rank " << ranks[i] );
    assertion ( ranks[i] >= 0 );
    ranksArray[i] = ranks[i];
  }
  // Create subgroup, containing processes contained in ranks
  preciceDebug ( "Restrict Group" );
  MPI_Group restrictedGroup;
  MPI_Group_incl ( currentGroup, ranks.size(), ranksArray, &restrictedGroup );
# ifdef Asserts
  int restrictedGroupSize = 0;
  MPI_Group_size ( restrictedGroup, & restrictedGroupSize );
  assertion ( restrictedGroupSize > 0 );
# endif
  // Create communicator, containing process of restrictedGroup
  preciceDebug ( "Create Comm" );
//  MPI_Comm restrictedCommunicator;
  MPI_Comm_create ( _globalCommunicator, restrictedGroup, &restrictedCommunicator );
  preciceDebug ( "Barrier" );
  MPI_Barrier ( _globalCommunicator );
  preciceDebug ( "Free current group" );
  MPI_Group_free ( &currentGroup );
  preciceDebug ( "Free restricted group" );
  MPI_Group_free ( &restrictedGroup );
  delete[] ranksArray;
# endif // not PRECICE_NO_MPI
  return restrictedCommunicator;
}

const std::vector<Parallel::AccessorGroup>& Parallel:: getAccessorGroups()
{
  preciceTrace ( "getAccessorGroups()" );
  assertion ( _isInitialized );
  return _accessorGroups;
}

}} // precice, utils

//#endif // not PRECICE_NO_MPI
