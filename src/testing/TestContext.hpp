#pragma once

#include <set>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "m2n/SharedPointer.hpp"
#include "precice/impl/Types.hpp"
#include "utils/Parallel.hpp"

namespace precice::testing {

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
  std::string_view name;

  /// the amount of ranks this participant runs on
  int size = 1;

  /// whether to initialize an intra-participant communication for this participant
  bool initIntraComm = false;

  /// Constructs a serial participant with a given name
  constexpr explicit ParticipantState(std::string_view n)
      : name(std::move(n)) {};

  /** Injects the amount of ranks this participant should run on.
   *
   * This call operator allows to write `"Fluid"_on(3_ranks)`
   *
   * @param[in] rsize the amount of Ranks to run on.
   *
   * @returns A reference to the ParticipantState allowing for chaining.
   */
  constexpr ParticipantState &operator()(Ranks rsize)
  {
    size = rsize.value;
    return *this;
  }

  /** Marks that this ParticipantState should initialize an intra-participant connection.
   *
   * @returns A reference to the ParticipantState allowing for chaining.
   */
  constexpr ParticipantState &setupIntraComm()
  {
    initIntraComm = true;
    return *this;
  }
};

/// User-defined literal allowing to create a serial ParticipantState from a given string.
inline constexpr ParticipantState operator""_on(const char *name, std::size_t)
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
  /// Ginkgo initialization
  Ginkgo,
};

/// Contains the setup description of a test including participants and requirements
struct TestSetup {

  /** Create a context representing an unnamed Participant running on a given count of Ranks and some requirements
   *
   * @note You need to construct a Participant if you require initializing
   * an intra-participant connection `"Serial"_on(3_ranks).setupIntraComm()`
   *
   * @see Require
   */
  template <class... T>
  TestSetup(T... args)
  {
    (handleOption(args), ...);
  }

  /// @{
  /// @name Option Handling
  void handleOption(ParticipantState participants);
  void handleOption(Ranks ranks);
  void handleOption(testing::Require requirement);
  /// @}

  /// total amount of ranks required by this setup
  int totalRanks() const;

  /// whether to initialize PETSc
  bool petsc = false;

  /// whether to initialize events
  bool events = false;

  /// whether to initialize Ginkgo (the device)
  bool ginkgo = false;

  /// All known participants
  std::vector<ParticipantState> participants;
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
 * @note Do not use this type directly. Use @ref PRECICE_TEST() and @ref PRECICE_TEST_SETUP() instead.
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

  /** Creates a context for a rank in the given TestSetup
   *
   * Unneeded ranks are marked as invalid.
   * Provides a TestContext named `context` which can be used in the test.
   *
   * @attention This call synchronizes all ranks
   *
   * @see @ref PRECICE_TEST()
   */
  TestContext(TestSetup setup);

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
   * @param[in] connector the requesting participant
   * @param[in] options a set of options concerning the created connection
   *
   * @note This function throws if the acceptor or connector are unknown!
   *
   * @see ConnectionOptions
   */
  m2n::PtrM2N connectPrimaryRanks(const std::string &acceptor, const std::string &connector, const ConnectionOptions &options = ConnectionOptions{}) const;

  /// Provides a user- and log-friendly description of the current context
  std::string describe() const;

private:
  TestSetup _setup;

  /// whether this context needs to initialize the intracomm
  bool _initIntraComm = false;

  /// the MPI communicator of the context
  utils::Parallel::CommStatePtr _contextComm;

  /// contains the name of every known Participant
  std::set<std::string> _names;

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

} // namespace precice::testing
