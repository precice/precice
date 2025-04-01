#include "Parallel.hpp"
#include <algorithm>
#include <map>
#include <memory>
#include <numeric>
#include <ostream>
#include <utility>
#include <vector>
#include "assertion.hpp"
#include "com/MPIDirectCommunication.hpp"
#include "logging/LogMacros.hpp"
#include "logging/Logger.hpp"

namespace precice::utils {

logging::Logger Parallel::_log("utils::Parallel");

Parallel::CommStatePtr Parallel::_currentState            = nullptr;
bool                   Parallel::_mpiInitializedByPrecice = false;

Parallel::InitializationState Parallel::_initState = Parallel::InitializationState::Uninitialized;

/// BEGIN CommState

Parallel::CommState::CommState(Parallel::CommState &&other) noexcept
{
  comm       = other.comm;
  other.comm = MPI_COMM_NULL;
  _owning    = other._owning;
  parent     = std::move(other.parent);
}

Parallel::CommState &Parallel::CommState::operator=(Parallel::CommState &&other) noexcept
{
  comm       = other.comm;
  other.comm = MPI_COMM_NULL;
  _owning    = other._owning;
  parent     = std::move(other.parent);
  return *this;
}

Parallel::CommState::~CommState() noexcept
{
#ifndef PRECICE_NO_MPI
  if (_owning && comm != MPI_COMM_SELF && comm != MPI_COMM_NULL && comm != MPI_COMM_WORLD) {
    MPI_Comm_free(&comm);
  }
#endif // not PRECICE_NO_MPI
}

int Parallel::CommState::rank() const
{
  int processRank = 0;
#ifndef PRECICE_NO_MPI
  MPI_Comm_rank(comm, &processRank);
#endif // not PRECICE_NO_MPI
  return processRank;
}

int Parallel::CommState::size() const
{

  int communicatorSize = 1;
#ifndef PRECICE_NO_MPI
  MPI_Comm_size(comm, &communicatorSize);
#endif // not PRECICE_NO_MPI
  return communicatorSize;
}

void Parallel::CommState::synchronize() const
{
#ifndef PRECICE_NO_MPI
  PRECICE_TRACE();
  if (!isNull()) {
    MPI_Barrier(comm);
  }
#endif // not PRECICE_NO_MPI
}

bool Parallel::CommState::isNull() const
{
#ifndef PRECICE_NO_MPI
  return comm == MPI_COMM_NULL;
#else
  return true;
#endif
}

Parallel::CommStatePtr Parallel::CommState::world()
{
#ifndef PRECICE_NO_MPI
  return fromComm(MPI_COMM_WORLD);
#else
  return null();
#endif
}

Parallel::CommStatePtr Parallel::CommState::null()
{
  return std::make_shared<CommState>();
}

Parallel::CommStatePtr Parallel::CommState::self()
{
#ifndef PRECICE_NO_MPI
  return fromComm(MPI_COMM_SELF);
#else
  return null();
#endif
}

Parallel::CommStatePtr Parallel::CommState::fromComm(Communicator comm)
{
  CommStatePtr state = null();
  state->comm        = comm;
  return state;
}

Parallel::CommStatePtr Parallel::CommState::fromExtern(Communicator comm)
{
  CommStatePtr state = null();
  state->comm        = comm;
  state->_owning     = false;
  return state;
}

void Parallel::CommState::print(std::ostream &out) const
{
  if (comm == MPI_COMM_NULL) {
    out << "COMM_NULL:invalid";
    return;
  }
#ifndef PRECICE_NO_MPI
  if (comm == MPI_COMM_SELF) {
    out << "COMM_SELF:1/1";
    return;
  }
  out << "COMM" << ((comm == MPI_COMM_WORLD) ? "_WORLD:" : ":");
  out << rank() << '/' << size();
  if (!_owning)
    out << "EXTERN";
#endif
}

/// END CommState

Parallel::CommStatePtr Parallel::current()
{
  if (!_currentState) {
    _currentState = CommState::world();
  }

  return _currentState;
}

void Parallel::resetCommState()
{
  _currentState = CommState::world();
}

void Parallel::resetManagedMPI()
{
  _mpiInitializedByPrecice = false;
  _initState               = Parallel::InitializationState::Uninitialized;
}

void Parallel::pushState(CommStatePtr newState)
{
  PRECICE_TRACE();
  PRECICE_ASSERT(newState != nullptr, "pushState cannot to be called with nullptr!");
  PRECICE_ASSERT(newState->parent == nullptr, "The parent of the given state must be empty!");
#ifndef NDEBUG
  PRECICE_DEBUG("Update comm state from {} to {}", *current(), *newState);
#endif
  newState->parent = _currentState;
  _currentState    = std::move(newState);
}

bool Parallel::isMPIInitialized()
{
#ifndef PRECICE_NO_MPI
  int isMPIInitialized{-1};
  MPI_Initialized(&isMPIInitialized);
  return isMPIInitialized != 0;
#else
  return false;
#endif // not PRECICE_NO_MPI
}

void Parallel::initializeOrDetectMPI(std::optional<Communicator> userProvided)
{
#ifndef PRECICE_NO_MPI
  PRECICE_ASSERT(!_mpiInitializedByPrecice,
                 "MPI cannot be initialized twice. You need to handle the MPI lifetime yourself.");

  bool isInit = isMPIInitialized();

  // Handle user-provided Communicator first
  if (userProvided.has_value()) {
    PRECICE_ASSERT(isInit, "A user-provided comm can only exist if MPI has been initialized.");
    pushState(Parallel::CommState::fromExtern(*userProvided));
    _initState = InitializationState::Provided;
    return;
  }

  // preCICE needs to initialize MPI itself
  if (!isInit) {
    MPI_Init(nullptr, nullptr);
    _currentState            = CommState::world();
    _initState               = InitializationState::Managed;
    _mpiInitializedByPrecice = true;
    return;
  }

  // preCICE is constructed in testing mode
  if (_currentState == nullptr) {
    // User initialized MPI but didn't pass a communicator
    // We need to use MPI_COMM_WORLD
    _currentState = CommState::world();
    _initState    = InitializationState::Unmanaged;
    return;
  }

  // We are in testing mode as \ref _currentState has been altered.
  _initState = InitializationState::Testing;
  return;
#endif
}

void Parallel::finalizeOrCleanupMPI()
{
#ifndef PRECICE_NO_MPI
  // Make sure all com states are freed at this point in time
  if (_initState == InitializationState::Testing) {
    resetCommState();
  } else {
    _currentState = nullptr;
  }

  if (_initState == InitializationState::Managed) {
    MPI_Finalize();
    PRECICE_ASSERT(_mpiInitializedByPrecice, "Something changed this state!");
  }

  _initState = InitializationState::Uninitialized;
#endif // not PRECICE_NO_MPI
}

void Parallel::initializeTestingMPI(
    int    *argc,
    char ***argv)
{
#ifndef PRECICE_NO_MPI
  PRECICE_ASSERT(!isMPIInitialized(), "MPI was already initialized.");
  MPI_Init(argc, argv);
  // By altering the commstate, preCICE will know that it is testing mode
  _currentState = CommState::world();
#endif // not PRECICE_NO_MPI
}

void Parallel::finalizeTestingMPI()
{
  // Make sure all com states are freed at this point in time
  resetCommState();
#ifndef PRECICE_NO_MPI
  MPI_Finalize();
#endif // not PRECICE_NO_MPI
}

// State altering

void Parallel::splitCommunicator(std::optional<int> group)
{
#ifndef PRECICE_NO_MPI
  PRECICE_TRACE(group.value_or(-1));

  // Passing the same key maintains rank order
  constexpr int keepTheSameOrder{0};

  MPI_Comm newComm = MPI_COMM_NULL;
  // Passing MPI_UNDEFINED as group results in MPI_COMM_NULL
  auto err = MPI_Comm_split(current()->comm, group.value_or(MPI_UNDEFINED), keepTheSameOrder, &newComm);
  PRECICE_ASSERT(err == MPI_SUCCESS, "MPI_Comm_split failed!", group.has_value(), group.value_or(-1));

  // Assemble and set new state
  pushState(CommState::fromComm(newComm));
#endif // not PRECICE_NO_MPI
}

void Parallel::popState()
{
  PRECICE_TRACE();
  auto state = current();
  if (state->parent) {
    _currentState = state->parent;
  }
}

int Parallel::getProcessRank()
{
  // Do not use TRACE or DEBUG here!
#ifndef PRECICE_NO_MPI
  if (_currentState) {
    return _currentState->rank();
  }
#endif // not PRECICE_NO_MPI
  return 0;
}

std::ostream &operator<<(std::ostream &out, const Parallel::CommState &value)
{
  value.print(out);
  return out;
}

} // namespace precice::utils

// #endif // not PRECICE_NO_MPI
