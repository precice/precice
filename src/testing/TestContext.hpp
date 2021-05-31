#pragma once

#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>
#include "m2n/M2N.hpp"
#include "m2n/SharedPointer.hpp"
#include "utils/Parallel.hpp"

namespace precice {
namespace testing {

/** Represents a count of MPI Ranks
 * @see TestContext()
 */
struct Ranks {
  int value;
};

/** User-defined literal for expressively defining multiple ranks
 *
 * @param[in] value the amount of ranks <= 1
 *
 * @returns a strong typed count of ranks
 */
inline constexpr Ranks operator""_ranks(unsigned long long value)
{
  return (value <= 1) ? throw std::runtime_error{"Cannot create a single rank with _ranks()! Use _rank() instead!"} : Ranks{static_cast<int>(value)};
}

/** User-defined literal for expressively defining a single rank
 *
 * @param[in] value the amount of ranks which has to be 1
 *
 * @returns a strong typed count of 1 rank
 */
inline constexpr Ranks operator""_rank(unsigned long long value)
{
  return (value == 1) ? Ranks{1} : throw std::runtime_error{"Cannot create multiple ranks with _rank()! Use _ranks() instead!"};
}

/// Represents a Partiticipant in a test
struct Participant {
  /// the name of the participant
  std::string name;

  /// the amount of ranks this participant runs on
  int size = 1;

  /// wheather to initialize a master-slave communication for this participant
  bool initMS = false;

  /// Constructs a serial participant with a given name
  explicit Participant(std::string n)
      : name(std::move(n)){};

  /** Injects the amount of ranks this participant should run on.
   *
   * This call operator allows to write `"Fluid"_on(3_ranks)`
   *
   * @param[in] rsize the amount of Ranks to run on.
   *
   * @returns A reference to the Participant allowing for chaining.
   */
  Participant &operator()(Ranks rsize)
  {
    size = rsize.value;
    return *this;
  }

  /** Marks that this Partiticipant should initialize a master-slave connection.
   *
   * @returns A reference to the Participant allowing for chaining.
   */
  Participant &setupMasterSlaves()
  {
    initMS = true;
    return *this;
  }
};

/// User-defined literal allowing to create a serial Participant from a given string.
inline Participant operator""_on(const char *name, std::size_t)
{
  return Participant{name};
}

static_assert(std::is_same<Participant &, decltype(""_on(1_rank))>::value, "");
static_assert(std::is_same<Participant &, decltype(""_on(2_ranks))>::value, "");

/** Defines requirements for a test setup
 *
 * @note These are used for unit-tests.
 * Integration tests calling the SolverInterface initialize required components themselves.
 */
enum class Require {
  /// Require to initialize PETSc. This implies the initialization of Events
  PETSc,
  /// Require to initialize Event.
  Events,
};

/** A type of distributed connection
 *
 * @see ConnectionOptions
 */
enum struct ConnectionType {
  GatherScatter,
  PointToPoint
};

/** Type representing options for an inter-participant connection.
 *
 * @see TestContext::connectMasters()
 * @see M2N::M2N()
 */
struct ConnectionOptions {
  ConnectionOptions() = default;

  /** Wheather to use only the Master-Master connection
   * @see M2N::M2N()
   */
  bool useOnlyMasterCom = false;

  /** Wheather to enable the two-level initialization
   * @see M2N::M2N()
   */
  bool useTwoLevelInit = false;

  /** The type of \ref DistributedCommunication to create
   * @see M2N::M2N()Q
   */
  ConnectionType type = ConnectionType::GatherScatter;
};

/** Type representing the context of a test.
 *
 * @note Do not use this type directly. Use @ref PRECICE_TEST() instead.
 *
 * This type is responsible for
 * 1. making sure that there are enough ranks to run the test on.
 * 2. restricting and splitting the MPI Communicator.
 * 3. handling invalid contexts (such as unneeded ranks)
 * 4. initializing the master-slave communication if requested. initializeMasterSlave()
 * 5. handling further requirements @see Require
 * 6. providing a usable context during the test isNamed(), isMaster(), isRank()
 * 7. cleaning up after the test case
 */
class TestContext {
public:
  using Participants = std::vector<Participant>;

  /// the name of the current participant
  std::string name;

  /// the rank of the current participant
  int rank = 0;

  /// the size of the Communicator of the current participant
  int size = 1;

  /// wheather this context is valid or not
  bool invalid = false;

  /// @{
  /// @name Construction

  /// Create a context representing an unnamed serial Participant
  TestContext() = default;

  /** Create a context representing an unnamed Participant running on a given count of Ranks
   *
   * @note You need to construct a Participant if you require initializing
   * a master-slave connection `"Serial"_on(3_ranks).setupMasterSlaves()`
   *
   * @attention This call synchonizes all ranks
   *
   */
  template <class... T>
  TestContext(Ranks ranks)
      : _simple(true)
  {
    Participants participants{"Serial"_on(ranks)};
    initialize(participants);
  }

  /** Create a context representing an unnamed Participant running on a given count of Ranks and some requirements
   *
   * @note You need to construct a Participant if you require initializing
   * a master-slave connection `"Serial"_on(3_ranks).setupMasterSlaves()`
   *
   * @attention This call synchonizes all ranks
   *
   * @see Require
   */
  template <class... T>
  TestContext(Ranks ranks, T... args)
      : _simple(true)
  {
    Participants participants{"Serial"_on(ranks)};
    handleOptions(participants, args...);
    initialize(participants);
  }

  /** Create a context representing one or more participants
   *
   * @attention This call synchonizes all ranks
   *
   * @see Require
   */
  template <class... T>
  TestContext(T... args)
  {
    Participants participants;
    handleOptions(participants, args...);
    initialize(participants);
  }

  /// @}

  /** Cleans-up all initialized parts and synchonizes all ranks
   * @attention This call synchonizes all ranks
   */
  ~TestContext() noexcept;

  /// Check wheater this context has a given size
  bool hasSize(int size) const;

  /// Check wheater this context has a given name
  bool isNamed(const std::string &name) const;

  /// Check wheater this context has a given rank inside the Partiticipant
  bool isRank(int rank) const;

  /** Check wheater this context is the master of a Participants
   * @note This is equivalent to `isRank(0)`
   */
  bool isMaster() const;

  /** Creates a M2N and establishes a master-master connection between participants
   * @param[in] acceptor the accepting participant
   * @param[in] requestor the requesting participant
   * @param[in] options a set of options concerning the created connection
   *
   * @note This function throws if the acceptor or requestor are unknown!
   *
   * @see ConnectionOptions
   */
  m2n::PtrM2N connectMasters(const std::string &acceptor, const std::string &requestor, const ConnectionOptions &options = ConnectionOptions{}) const;

  /// Provides a user- and log-friendly description of the current context
  std::string describe() const;

private:
  /// wheater to initialize PETSc
  bool _petsc = false;

  /// wheater to initialize events
  bool _events = false;

  /// wheater this Context was created with a Ranks constructor
  bool _simple = false;

  /// wheater to initialize a master-slave connection
  bool _initMS = false;

  /// the MPI communicator of the context
  utils::Parallel::CommStatePtr _contextComm;

  /// contains the name of every known Participant
  std::vector<std::string> _names;

  /// @{
  /// @name Option Handling
  void handleOption(Participants &participants, Participant participant);
  void handleOption(Participants &participants, testing::Require requirement);

  template <class LastOption>
  void handleOptions(Participants &participants, LastOption &last)
  {
    handleOption(participants, last);
  }

  template <class NextOption, class... Rest>
  void handleOptions(Participants &participants, NextOption &next, Rest &... rest)
  {
    handleOption(participants, next);
    handleOptions(participants, rest...);
  }
  /// @}

  /** set the context from a Participants and a given rank
   * Both uniquely identify a context.
   */
  void setContextFrom(const Participant &p, int rank);

  /// @{
  /// @name Initialization

  /// Main entrypoint
  void initialize(const Participants &participants);

  /** Check, restrict and split the MPI communicator
   * Marks unneeded contexts as invalid
   */
  void initializeMPI(const Participants &participants);

  /// Initialize the Master-Slave connection if requested
  void initializeMasterSlave();

  /// Initialize PETSc if required
  void initializePetsc();

  /// Initialize Events if required
  void initializeEvents();

  /// @}
};

} // namespace testing
} // namespace precice
