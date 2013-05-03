// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_UTILS_HELPERSPARALLEL_HPP_
#define PRECICE_UTILS_HELPERSPARALLEL_HPP_

#include "tarch/logging/Log.h"
#include <vector>
#include <string>

#ifndef PRECICE_NO_MPI

#include "mpi.h"

#define PRECICE_MASTER_ONLY \
   if ( precice::utils::Parallel::getProcessRank() == 0 )
#else
#define PRECICE_MASTER_ONLY
#endif // not PRECICE_NO_MPI


//#ifndef PRECICE_NO_MPI

namespace precice {
namespace utils {

/**
 * @brief Utility class for managing MPI operations.
 */
class Parallel
{
public:

  #ifndef PRECICE_NO_MPI
  typedef MPI_Comm Communicator;
  #else
  typedef int Communicator;
  #define MPI_COMM_NULL -1
  #endif

  /**
   * @brief Used to sort and order all coupling participants.
   */
//  struct Accessor {
//    std::string name;
//    int groupID;
//    int leaderRank;
//    int groupSize;
//  };

  struct AccessorGroup {
    std::string name;
    int id;
    int leaderRank;
    int size;
  };

  static Communicator getCommunicatorWorld();

  static void initialize (
    int*               argc,
    char***            argv,
    const std::string& groupName );

  /**
   * @brief Finalizes MPI environment.
   */
  static void finalize();

  /**
   * @brief Returns the global process rank.
   */
  static int getProcessRank();

  /**
   * @brief Returns the local process rank.
   *
   * If only one accessor group is present, returns getProcessRank().
   */
  static int getLocalProcessRank();

  /**
   * @brief Returns the number of processes in the communicator.
   */
  static int getCommunicatorSize();

  /**
   * @brief Synchronizes all processes.
   */
  static void synchronizeProcesses();

  /**
   * @brief Synchronizes all local processes.
   *
   * If only one accessor group is present, calls synchcronizeProcesses().
   */
  static void synchronizeLocalProcesses();

  /**
   * @brief Switches precice communication away from global space to given one.
   *
   * The switch has only effects on communication means created after the
   * switch. The ones before stay in their old communication universe.
   * Standard communication space is MPI_COMM_WORLD. The local process rank
   * and communicator size is recomputed, relative to the new default
   * communicator.
   *
   * ATTENTION: Will result in an error, if called by a process not in the new
   *            default communicator!
   */
  static void setGlobalCommunicator ( Communicator defaultCommunicator );

  /**
   * @brief Returns the default communicator.
   */
  static const Communicator& getGlobalCommunicator();

  /**
   * @brief Returns communicator of processes within one group.
   *
   * This communicator is empty, until initialize() has been called.
   */
  static const Communicator& getLocalCommunicator();

  /**
   * @brief Returns a communicator with a subset of processes.
   *
   * Does not change the default communicator.
   *
   * ATTENTION: Has to be called by every process in the communicator to be
   *            restricted, otherwise, a deadlock is achieved!
   *
   * @param ids [IN] Process ranks to be selected for restricted comm.
   */
  static Communicator getRestrictedCommunicator ( const std::vector<int>& ranks );

//  static void exchangeGroupInformation ( const std::string & groupName );

  static const std::vector<AccessorGroup>& getAccessorGroups();

private:

  static tarch::logging::Log _log;

  static Communicator _globalCommunicator;

  static Communicator _localCommunicator;

  //static bool _isLocalCommunicatorSet;

  // @brief Processes participating in direct communication.
  static std::vector<AccessorGroup> _accessorGroups;

  static bool _isInitialized;

  static bool _autoInitialized;
};


}} // namespace precice, utils

//#endif // not PRECICE_NO_MPI

#endif /* PRECICE_UTILS_HELPERSPARALLEL_HPP_ */
