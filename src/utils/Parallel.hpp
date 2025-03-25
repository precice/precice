#pragma once

#include <iosfwd>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "logging/Logger.hpp"
#include "precice/impl/Types.hpp"

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

  struct CommState;

  using CommStatePtr = std::shared_ptr<CommState>;

  /** Represents a Communicator state based on a parent state
   * An object of this type owns it communicator and will free it at the end of its lifetime.
   * It also owns the CommState it originated from (parent) via shared ownership.
   */
  struct CommState {
    /// The native communicator that represents this state
    Communicator comm = MPI_COMM_NULL;

    /// A shared pointer to the parent CommState
    CommStatePtr parent = nullptr;

    /// Whether this state owns the communicator and has to free it.
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
    Rank rank() const;

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

  /// Return true if MPI is initialized
  static bool isMPIInitialized();

  /** Initializes or detects an existing MPI environment
   *
   * If a custom MPI Communicator is provided via \ref userProvided then this registers a user-provided MPI session.
   *
   * If MPI has already been initialized, then preCICE uses MPI_COMM_WORLD as communicator
   * and registers a unmanaged MPI session.
   * If \ref _currentState isn't nullptr, then this signals launch inside a test.
   *
   * If MPI hasn't been initialized yet, then preCICE takes ownership.
   * It initializes the environment and will later destroy it.
   * As MPI forbids reinitialization, this prevents reconstruction.
   *
   * @param[in] userProvided an optional user-provided Communicator
   *
   * @see finalizeOrCleanupMPI()
   */
  static void initializeOrDetectMPI(std::optional<Communicator> userProvided = std::nullopt);

  /**
   * @brief Finalized a managed MPI environment or cleans up after an non-managed session.
   *
   * To initialize the managed MPI call @ref initializeOrDetectMPI().
   * This finalizes MPI only if the preCICE initialized the MPI session itself.
   *
   * @see initializeOrDetectMPI()
   */
  static void finalizeOrCleanupMPI();

  /** Unconditionally initializes the MPI environment.
   *
   * Alters the \ref _currentState, which indicates a testing session.
   *
   * @param[in] argc Parameter count
   * @param[in] argv Parameter values, is passed to MPI_Init
   */
  static void initializeTestingMPI(
      int *   argc,
      char ***argv);

  /// Unconditionally finalizes MPI environment.
  static void finalizeTestingMPI();

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
  static void splitCommunicator(std::optional<int> group = std::nullopt);

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

  /// @name Misc
  /// @{

  /** Returns the global process rank or 0
   * used in assertions.
   */
  static Rank getProcessRank();

  /// Returns an owning pointer to the current CommState.
  static CommStatePtr current();

  /// @}

private:
  static logging::Logger _log;

  static CommStatePtr _currentState;

  /// Flag to safeguard against reinitializing MPI, which is forbidden
  static bool _mpiInitializedByPrecice;

  /// Kind of initialization that took place
  enum struct InitializationState {
    Uninitialized, /// Not initialized
    Provided,      /// Communicator was provided by the user
    Managed,       /// preCICE manages the lifetime of the MPI environment
    Unmanaged,     /// preCICE was initialized in an existing MPI environment
    Testing        /// preCICE was initialized in a testing environment initialized with \ref initializeTestingMPI()
  };

  static InitializationState _initState;

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

template <>
struct fmt::formatter<precice::utils::Parallel::CommState> : ostream_formatter {
};
