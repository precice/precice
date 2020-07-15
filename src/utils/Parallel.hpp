#pragma once

#include <iosfwd>
#include <memory>
#include <string>
#include <vector>
#include "logging/Logger.hpp"

#ifndef PRECICE_NO_MPI
#include <mpi.h>
#endif // not PRECICE_NO_MPI

namespace precice {
namespace logging {
class Logger;
} // namespace logging

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

  struct CommState;

  using CommStatePtr = std::shared_ptr<CommState>;

  /** Represents a Communicator state based on a parent state
   * An object of this type owns it communicator and will free it at the end of its lifetime.
   * It also owns the CommState it originated from (parent) via shared ownership.
   */
  struct CommState {

    /// The different groups that were used to split the communicator
    std::vector<AccessorGroup> groups;

    /// The native communicator that represents this state
    Communicator comm = MPI_COMM_NULL;

    /// A shared pointer to the parent CommState
    CommStatePtr parent = nullptr;

    /// Wether this state owns the communicator and has to free it.
    bool _owning = true;

    /// @name Construction and Destruction
    /// @{
    CommState() = default;

    CommState(const CommState &) = delete;
    CommState &operator=(const CommState &) = delete;

    CommState(CommState &&) noexcept;
    CommState &operator=(CommState &&) noexcept;

    /// Frees the communicator if allowed
    ~CommState() noexcept;

    /// @}

    /// @name Factory functions
    /// @{

    /// returns a commstate containing MPI_COMM_WORLD
    static CommStatePtr world();

    /// returns an blank commstate representing MPI_COMM_NULL
    static CommStatePtr null();

    /// returns the commstate representing MPI_COMM_SELF
    static CommStatePtr self();

    /// returns the commstate representing comm
    static CommStatePtr fromComm(Communicator comm);

    /** returns the commstate representing an extern comm
     *
     * @attention This state is not owning and does not call free on comm!
     *
     * @see _owning
     */
    static CommStatePtr fromExtern(Communicator comm);

    /// @}

    /// @name Communicator Access
    /// @{

    /// Returns the current rank in comm
    int rank() const;

    /// Returns size of comm
    int size() const;

    /// Returns weather the comm is NULL
    bool isNull() const;

    /// @}

    /// @name Misc
    /// @{

    /** Synchronizes all processes in the communicator
     * @attention This is a collective operation and has to be called by every rank in the communicator comm!
     */
    void synchronize() const;

    /// pretty printer for comms
    void print(std::ostream &out) const;

    /// @}
  };

  /// @name Initialization and Finalization
  /// @{

  /**
   * @brief Initializes the MPI environment and manages it.
   *
   * This keeps track of the state when called by setting _isInitialized and _mpiInitializedByPrecice
   * To finalize the managed MPI call @ref finalizeManagedMPI().
   *
   * @param[in] argc Parameter count
   * @param[in] argv Parameter values, is passed to MPI_Init
   *
   * @see finalizeManagedMPI
   */
  static void initializeManagedMPI(
      int *   argc,
      char ***argv);

  /**
   * @brief Unconditionally initializes the MPI environment.
   *
   * @param[in] argc Parameter count
   * @param[in] argv Parameter values, is passed to MPI_Init
   */
  static void initializeMPI(
      int *   argc,
      char ***argv);

  /**
   * @brief Finalized a managed MPI environment.
   *
   * To initialize the managed MPI call @ref initializeManagedMPI().
   * This finalizes MPI only if it was not initialized before the call to @ref initializeManagedMPI()
   *
   * @see InitializeManagedMPI
   */
  static void finalizeManagedMPI();

  /// Unconditionally finalizes MPI environment.
  static void finalizeMPI();

  /// Registers a user-provided communicator
  static void registerUserProvidedComm(Communicator comm);

  /// @}

  /// @name State-altering Functions
  /// @{

  /**
   * @brief Splits and creates a local MPI communicator according to groupName
   *
   * Updates the current CommState of Parallel.
   *
   * @param[in] groupName MPI group in which the calling process will be put in
   */
  static void splitCommunicator(const std::string &groupName);

  /// Creates a restricted communicator
  /**
   * Set the new, restricted communicator on all ranks that are contained in
   * that new communicator. Sets the communicator of the other ranks to MPI_COMM_NULL.
   *
   * @param[in] newSize How many ranks to restrict the Communicator to
   *
   * @postcondition current()->size() == newSize || current()->isNull()
   */
  static void restrictCommunicator(int newSize);

  /** Resets the commState to World
   *
   * This resets the _currentState to CommState::world().
   * Unrequired Communicators will be automatically freed.
   *
   * @postcondition cureent() == CommState::world()
   * 
   */
  static void resetCommState();

  /** Resets managed MPI
   *
   * This resets _mpiInitializedByPrecice and _isInitialized to false
   */
  static void resetManagedMPI();

  /** Sets the current state to its parent
   *
   * @postcondition current() == old current()->parent
   */
  static void popState();
  /// @}

  /// @name State Access
  /// @{

  /// clears groups for communicator splitting
  // @todo remove
  // static void clearGroups(){};

  /// Returns the global process rank.
  //@todo remove
  static int getProcessRank();

  /**
   * @brief Returns the local process rank.
   *
   * If only one accessor group is present, returns getProcessRank().
   */
  //@todo remove
  static int getLocalProcessRank();

  /// Returns the number of processes in the global communicator.
  //@todo remove
  // static int getCommunicatorSize();

  /// Returns the number of processes in the given communicator.
  //@todo remove
  // static int getCommunicatorSize(Communicator comm);

  /// Synchronizes all processes.
  //@todo remove
  // static void synchronizeProcesses();

  /**
   * @brief Synchronizes all local processes.
   *
   * If only one accessor group is present, calls synchcronizeProcesses().
   */
  //@todo remove
  // static void synchronizeLocalProcesses();

  /**
   * @brief Switches precice communication away from global space to given one.
   *
   * The switch has only effects on communication means created after the
   * switch. The ones before stay in their old communication universe.
   * Standard communication space is MPI_COMM_WORLD. The local process rank
   * and communicator size is recomputed, relative to the new default
   * communicator.
   *
   * @param[in] defaultCommunicator The new global/default Communicator
   * @param[in] free free the old communicator?
   *
   * @attention Will result in an error, if called by a process not in the new
   *            default communicator!
   */
  //static void setGlobalCommunicator(Communicator defaultCommunicator, bool free = true);

  /// @}

  /// @name Misc
  /// @{

  /** Returns an owning pointer to the global CommState, being the parent of the current CommState
   *
   * @note Calling this on World returns World.
   * 
   * @see getLocalCommunicator()
   */
  static const CommStatePtr getGlobalCommState();

  /**
   * @brief Returns an owning pointer to the local CommState, being the current CommState
   *
   * @note equivalent to calling current()
   *
   * @see getGlobalCommunicator()
   */
  static const CommStatePtr getLocalCommState();

  /// Returns an owning pointer to the current CommState.
  static CommStatePtr current();

  /// @}

private:
  static logging::Logger _log;

  static CommStatePtr _currentState;

  static bool _isInitialized;

  static bool _mpiInitializedByPrecice;

  /** Pushes a new state on the state stack
   *
   * Sets the parent of newState to the current state.
   * Sets the current state to the newState.
   *
   * @precondition newState->parent == nullptr
   * @postcondition _currentState == newState
   * @postcondition _currentState->parent == old _currentState
   */
  static void pushState(CommStatePtr newState);
};

/// pretty printer for CommState
std::ostream &operator<<(std::ostream &out, const Parallel::CommState &value);

} // namespace utils
} // namespace precice
