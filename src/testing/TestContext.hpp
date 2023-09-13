#pragma once

#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "m2n/SharedPointer.hpp"
#include "precice/impl/Types.hpp"
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

/// Represents a ParticipantState in a test
struct ParticipantState {
  /// the name of the participant
  std::string name;

  /// the amount of ranks this participant runs on
  int size = 1;

  /// whether to initialize an intra-participant communication for this participant
  bool initIntraComm = false;

  /// Constructs a serial participant with a given name
  explicit ParticipantState(std::string n)
      : name(std::move(n)){};

  /** Injects the amount of ranks this participant should run on.
   *
   * This call operator allows to write `"Fluid"_on(3_ranks)`
   *
   * @param[in] rsize the amount of Ranks to run on.
   *
   * @returns A reference to the ParticipantState allowing for chaining.
   */
  ParticipantState &operator()(Ranks rsize)
  {
    size = rsize.value;
    return *this;
  }

  /** Marks that this ParticipantState should initialize an intra-participant connection.
   *
   * @returns A reference to the ParticipantState allowing for chaining.
   */
  ParticipantState &setupIntraComm()
  {
    initIntraComm = true;
    return *this;
  }
};

/// User-defined literal allowing to create a serial ParticipantState from a given string.
inline ParticipantState operator""_on(const char *name, std::size_t)
{
  return ParticipantState{name};
}

static_assert(std::is_same<ParticipantState &, decltype(""_on(1_rank))>::value, "");
static_assert(std::is_same<ParticipantState &, decltype(""_on(2_ranks))>::value, "");

/** Defines requirements for a test setup
 *
 * @note These are used for unit-tests.
 * Integration tests calling the Participant initialize required components themselves.
 */
enum class Require {
  /// Require to initialize PETSc. This implies the initialization of Events
  PETSc,
  /// Require to initialize Event.
  Events,
  /// Kokkos initialization
  Kokkos,
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
 * @see TestContext::connectPrimaryRanks()
 * @see M2N::M2N()
 */
struct ConnectionOptions {
  ConnectionOptions() = default;

  /** Whether to use only the primary connection
   * @see M2N::M2N()
   */
  bool useOnlyPrimaryCom = false;

  /** Whether to enable the two-level initialization
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
 * 4. initializing the intra-participant communication if requested. initializeIntraComm()
 * 5. handling further requirements @see Require
 * 6. providing a usable context during the test isNamed(), isPrimary(), isRank()
 * 7. cleaning up after the test case
 */
class TestContext {
public:
  using Participants = std::vector<ParticipantState>;

  /// the name of the current participant
  std::string name;

  /// the rank of the current participant
  Rank rank = 0;

  /// the size of the Communicator of the current participant
  int size = 1;

  /// whether this context is valid or not
  bool invalid = false;

  /// @{
  /// @name Construction

  /// Create a context representing an unnamed serial Participant
  TestContext() = default;

  /** Create a context representing an unnamed Participant running on a given count of Ranks
   *
   * @note You need to construct a Participant if you require initializing
   * an intra-participant connection `"Serial"_on(3_ranks).setupIntraComm()`
   *
   * @attention This call synchronizes all ranks
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
   * an intra-participant connection `"Serial"_on(3_ranks).setupIntraComm()`
   *
   * @attention This call synchronizes all ranks
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
   * @attention This call synchronizes all ranks
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

  /** Cleans-up all initialized parts and synchronizes all ranks
   * @attention This call synchronizes all ranks
   */
  ~TestContext() noexcept;

  /** Returns the canonical config name of this test.
   *
   * The location of integration tests are tied the test name and test suites.
   * This computes the canonical filename of this test's configuration file
   * based on the repository root, the current test suites and name.
   *
   * @return the full filepath of this test's configuration file
   */
  std::string config() const;

  /** Prefixes the given filename with the test directory.
   *
   * The filename will be located in the same directory as the current test file.
   *
   * @return the full filepath to the filename relative
   */
  std::string prefix(const std::string &filename) const;

  /// Check whether this context has a given size
  bool hasSize(int size) const;

  /// Check whether this context has a given name
  bool isNamed(const std::string &name) const;

  /// Check whether this context has a given rank inside the Participant
  bool isRank(Rank rank) const;

  /// Returns a pointer to the MPI communicator of this context
  auto comm()
  {
    return &(_contextComm->comm);
  }

  /** Check whether this context is the primary rank of a participant
   * @note This is equivalent to `isRank(0)`
   */
  bool isPrimary() const;

  /** Creates a M2N and establishes a primary connection between participants
   * @param[in] acceptor the accepting participant
   * @param[in] requestor the requesting participant
   * @param[in] options a set of options concerning the created connection
   *
   * @note This function throws if the acceptor or requestor are unknown!
   *
   * @see ConnectionOptions
   */
  m2n::PtrM2N connectPrimaryRanks(const std::string &acceptor, const std::string &requestor, const ConnectionOptions &options = ConnectionOptions{}) const;

  /// Provides a user- and log-friendly description of the current context
  std::string describe() const;

private:
  /// whether to initialize PETSc
  bool _petsc = false;

  /// whether to initialize events
  bool _events = false;

  /// whether to initialize PETSc
  bool _kokkos = false;

  /// whether this Context was created with a Ranks constructor
  bool _simple = false;

  /// whether to initialize an intra-participant connection
  bool _initIntraComm = false;

  /// the MPI communicator of the context
  utils::Parallel::CommStatePtr _contextComm;

  /// contains the name of every known Participant
  std::vector<std::string> _names;

  /// @{
  /// @name Option Handling
  void handleOption(Participants &participants, ParticipantState participant);
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

  /** set the context from a Participants and the current com
   * Both uniquely identify a context.
   * @see Par::current()
   */
  void setContextFrom(const ParticipantState &p);

  /// @{
  /// @name Initialization

  /// Main entrypoint
  void initialize(const Participants &participants);

  /** Check, restrict and split the MPI communicator
   * Marks unneeded contexts as invalid
   */
  void initializeMPI(const Participants &participants);

  /// Initialize the intra-participant communication connection if requested
  void initializeIntraComm();

  /// Initialize PETSc if required
  void initializePetsc();

  /// Initialize Events if required
  void initializeEvents();

  /// Initialize Ginkgo if required
  void initializeGinkgo();
  /// @}
};

} // namespace testing
} // namespace precice
