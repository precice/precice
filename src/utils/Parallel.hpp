#pragma once

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

namespace precice {
namespace utils {

/// Utility class for managing MPI operations.
class Parallel
{
public:

#ifndef PRECICE_NO_MPI
  typedef MPI_Comm Communicator;
#else
  typedef int Communicator;
#define MPI_COMM_NULL -1
#endif

  /// Used to sort and order all coupling participants.
  struct AccessorGroup {
    std::string name;
    int id;
    int leaderRank;
    int size;
  };

  static Communicator getCommunicatorWorld();

  /**
   * @brief Splits and creates a local MPI communicator according to groupName
   *
   * @param[in] groupName MPI group in which the calling process will be put in
   */
  static void splitCommunicator (
    const std::string& groupName );


  /**
   * @brief Initializes the MPI environment.
   *
   * @param[in] argc Parameter count
   * @param[in] argc Parameter values, is passed to MPI_Init
   */
  static void initializeMPI (
    int*               argc,
    char***            argv);

  /// Finalizes MPI environment.
  static void finalizeMPI();

  /// clears groups for communicator splitting
  static void clearGroups();

  /// Returns the global process rank.
  static int getProcessRank();

  /**
   * @brief Returns the local process rank.
   *
   * If only one accessor group is present, returns getProcessRank().
   */
  static int getLocalProcessRank();

  /// Returns the number of processes in the communicator.
  static int getCommunicatorSize();
  
  /// Synchronizes all processes.
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

  /// Returns the default communicator.
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

  static const std::vector<AccessorGroup>& getAccessorGroups();

private:

  static tarch::logging::Log _log;

  static Communicator _globalCommunicator;

  static Communicator _localCommunicator;

  /// Processes participating in direct communication.
  static std::vector<AccessorGroup> _accessorGroups;

  static bool _isInitialized;

  static bool _isSplit;
};


}} // namespace precice, utils
