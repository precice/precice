#pragma once

#include <string>
#include <vector>
#include "logging/Logger.hpp"

#ifndef PRECICE_NO_MPI
#include <mpi.h>
#endif // not PRECICE_NO_MPI

namespace precice {
namespace utils {

/// Utility class for managing MPI operations.
class Parallel {
public:
#ifndef PRECICE_NO_MPI
  using Communicator = MPI_Comm;
#else
  using Communicator = std::nullptr_t;
#define MPI_COMM_NULL nullptr
#endif

  /// Used to sort and order all coupling participants.
  struct AccessorGroup {
    std::string name;
    int         id;
    int         leaderRank;
    int         size;
  };

  static Communicator getCommunicatorWorld();

  /**
   * @brief Splits and creates a local MPI communicator according to groupName
   *
   * @param[in] groupName MPI group in which the calling process will be put in
   */
  static void splitCommunicator(const std::string &groupName);

  /**
   * @brief Initializes the MPI environment.
   *
   * @param[in] argc Parameter count
   * @param[in] argv Parameter values, is passed to MPI_Init
   */
  static void initializeMPI(
      int *   argc,
      char ***argv);

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
   * @attention Will result in an error, if called by a process not in the new
   *            default communicator!
   */
  static void setGlobalCommunicator(Communicator defaultCommunicator);

  /// Returns the default communicator.
  static const Communicator &getGlobalCommunicator();

  /**
   * @brief Returns communicator of processes within one group.
   *
   * This communicator is empty, until initialize() has been called.
   */
  static const Communicator &getLocalCommunicator();

  /**
   * @brief Returns a communicator with a subset of processes.
   *
   * Does not change the default communicator.
   *
   * @attention Has to be called by every process in the communicator to be
   *            restricted, otherwise, a deadlock is achieved!
   *
   * @param[in] ranks Process ranks to be selected for restricted comm.
   */
  static Communicator getRestrictedCommunicator(const std::vector<int> &ranks);

  /// Create a restricted communicator and sets them as the global communicator
  /**
   * Set the new, restricted communicator on all ranks that are contained in
   * that new communicator. Leaves the other ranks untouched.
   *
   * @param[in] ranks Process ranks to be selected for restricted comm.
   *
   */
  static void restrictGlobalCommunicator(const std::vector<int> &ranks);

  static const std::vector<AccessorGroup> &getAccessorGroups();

private:
  static logging::Logger _log;

  static Communicator _globalCommunicator;

  static Communicator _localCommunicator;

  /// Processes participating in direct communication.
  static std::vector<AccessorGroup> _accessorGroups;

  static bool _isInitialized;

  static bool _isSplit;

  static bool _mpiInitializedByPrecice;
};
} // namespace utils
} // namespace precice
