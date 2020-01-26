#include "Parallel.hpp"
#include <map>
#include "assertion.hpp"
#include "com/MPIDirectCommunication.hpp"

namespace precice {
namespace utils {

logging::Logger Parallel::_log("utils::Parallel");

bool                   Parallel::_isInitialized           = false;
Parallel::CommStatePtr Parallel::_currentState            = nullptr;
bool                   Parallel::_mpiInitializedByPrecice = false;

/// BEGIN CommState

Parallel::CommState::CommState(Parallel::CommState &&other) noexcept
{
  groups     = std::move(other.groups);
  comm       = other.comm;
  other.comm = MPI_COMM_NULL;
  parent     = std::move(other.parent);
}

Parallel::CommState &Parallel::CommState::operator=(Parallel::CommState &&other) noexcept
{
  groups     = std::move(other.groups);
  comm       = other.comm;
  other.comm = MPI_COMM_NULL;
  parent     = std::move(other.parent);
}

Parallel::CommState::~CommState() noexcept
{
  if (comm != MPI_COMM_SELF && comm != MPI_COMM_NULL && comm != MPI_COMM_WORLD) {
    MPI_Comm_free(&comm);
  }
}

int Parallel::CommState::rank() const
{
  int processRank = 0;
#ifndef PRECICE_NO_MPI
  if (_isInitialized) {
    MPI_Comm_rank(comm, &processRank);
  }
#endif // not PRECICE_NO_MPI
  return processRank;
}

int Parallel::CommState::size() const
{

  int communicatorSize = 1;
#ifndef PRECICE_NO_MPI
  if (_isInitialized) {
    MPI_Comm_size(comm, &communicatorSize);
  }
#endif // not PRECICE_NO_MPI
  return communicatorSize;
}

void Parallel::CommState::synchronize() const
{
#ifndef PRECICE_NO_MPI
  PRECICE_TRACE();
  PRECICE_ASSERT(_isInitialized);
  MPI_Barrier(comm);
#endif // not PRECICE_NO_MPI
}

bool Parallel::CommState::isNull() const
{
  return comm == MPI_COMM_NULL;
}

static Parallel::CommStatePtr Parallel::CommState::world()
{
  return fromComm(MPI_COMM_WORLD);
}

static Parallel::CommStatePtr Parallel::CommState::null()
{
  return std::make_shared<CommState>();
}

static Parallel::CommState::CommStatePtr Parallel::CommState::self()
{
  return fromComm(MPI_COMM_SELF);
}

static Parallel::CommState::CommStatePtr Parallel::CommState::fromComm(Communicator comm)
{
  CommStatePtr state = null();
  state->comm        = comm;
  return state;
}

std::ostream &Parallel::CommState::operator<<(std::ostream &out) const
{
  if (comm == MPI_COMM_SELF) {
    return out << "COMM_SELF:1/1";
  }
  if (comm == MPI_COMM_NULL) {
    return out << "COMM_NULL:invalid";
  }
  out << "COMM" << ((comm == MPI_COMM_WORLD) ? "_WORLD:" : ":");
  return out << rank() << '/' << size();
}

/// END CommState

static Parallel::CommStatePtr current()
{
  if (_currentState) {
    return _currentState;
  }

  _currentState = CommState::world();
}

static CommStatePtr resetCommState()
{
  _currentState = CommState::world();
}

static void Parallel::pushState(CommStatePtr newState)
{
  PRECICE_TRACE();
  PRECICE_ASSERT(newState, "pushState cannot to be called with nullptr!");
  PRECICE_ASSERT(newState->parent == false, "The parent of the given state must be empty!");
#ifndef NDEBUG
  PRECICE_DEBUG("Update comm state from " << *current() << " to " << *newState);
#endif
  newState->parent = _currentState;
  _currentState    = std::move(newState);
}

Parallel::Communicator Parallel::getCommunicatorWorld()
{
#ifndef PRECICE_NO_MPI
  return MPI_COMM_WORLD;
#else
  return nullptr;
#endif
}

void Parallel::initializeMPI(
    int *   argc,
    char ***argv)
{
#ifndef PRECICE_NO_MPI
  PRECICE_TRACE();
  int isMPIInitialized;
  MPI_Initialized(&isMPIInitialized);
  if (not isMPIInitialized) {
    PRECICE_DEBUG("Initialize MPI");
    _mpiInitializedByPrecice = true;
    MPI_Init(argc, argv);
  }
  _isInitialized = true;
#endif // not PRECICE_NO_MPI
}

void Parallel::splitCommunicator(const std::string &groupName)
{
#ifndef PRECICE_NO_MPI
  PRECICE_TRACE(groupName);

  // Exchange group information
  PRECICE_DEBUG("Exchange group information");
  //_accessorGroups.clear(); // Makes reinitialization possible

  CommStatePtr baseState  = current();
  MPI_Comm     globalComm = current->comm;
  int          rank       = baseState->rank();
  int          size       = baseState->size();

  PRECICE_ASSERT(size > 1, "Splitting a communicator or size 1 is not possible!");

  std::vector<AccessorGroup>  accessorGroups;
  com::MPIDirectCommunication com;

  if (rank == 0) { // Master
    // Step 1 receive groupNames from Salves and setup accessorGroups
    std::map<std::string, int> groupMap; // map from names to group ID
    groupMap[groupName] = 0;
    AccessorGroup newGroup;
    newGroup.id         = 0;
    newGroup.size       = 1;
    newGroup.leaderRank = 0;
    newGroup.name       = groupName;
    accessorGroups.push_back(newGroup);
    for (int i = 1; i < size; i++) {
      std::string name;
      com.receive(name, i); // Receive group name from all ranks
      if (groupMap.find(name) == groupMap.end()) {
        groupMap[name] = accessorGroups.size();
        AccessorGroup newGroup;
        newGroup.id         = accessorGroups.size();
        newGroup.size       = 1;
        newGroup.leaderRank = i;
        newGroup.name       = name;
        accessorGroups.push_back(newGroup);
      } else {
        accessorGroups[groupMap[name]].size++;
      }
    }
    // Step 2 send groupCount to Master
    auto groupCount = static_cast<int>(accessorGroups.size());
    MPI_Bcast(&groupCount, 1, MPI_INT, 0, globalComm);
    PRECICE_ASSERT(groupCount > 1, "Calling split with a single group is not permitted!");

    // Step 3 send AccessorGroups to Slaves
    for (const AccessorGroup &group : accessorGroups) {
      // @TODO can we use broadcast as the master sends this to everyone else?
      for (int i = 1; i < size; i++) {
        com.send(group.name, i);
        com.send(group.leaderRank, i);
        com.send(group.id, i);
        com.send(group.size, i);
      }
    }
  } else { // rank != 0
    // Step 1 send groupName to Master
    com.send(groupName, 0);
    // Step 2 receive groupCount from Master
    int groupCount = -1;
    MPI_Bcast(&groupCount, 1, MPI_INT, 0, globalComm);
    PRECICE_ASSERT(groupCount > 1, "Calling split with a single group is not permitted!");

    // Step 3 receive AccessorGroups from Master
    for (int i = 0; i < groupCount; i++) {
      AccessorGroup newGroup;
      com.receive(newGroup.name, 0);
      com.receive(newGroup.leaderRank, 0);
      com.receive(newGroup.id, 0);
      com.receive(newGroup.size, 0);
      accessorGroups.emplace_back(std::move(newGroup));
    }
  }

  // Step 4 split into groups
  PRECICE_DEBUG("Split groups");
  auto thisGroup = std::find_if(accessorGroups.begin(), accessorGroups.end(), [groupName](const AccessorGroup &group) { return group.name == groupName; });
  PRECICE_ASSERT(thisGroup != accessorGroups.end(), "This requested groupName \"" << groupName << "\" is not in accessorGroups!");
  // Create a new communicator that contains only ranks of my group
  MPI_Comm newComm = MPI_COMM_NULL;
  MPI_Comm_split(globalComm, thisGroup->id, rank, &newComm);
  // Assemble and set new state
  auto newState            = CommState::null();
  newState->comm           = newComm;
  newState->accessorGroups = std::move(accessorGroups);
  newState->parent         = _currentState;
  _currentState            = std::move(newState);

#ifndef NDEBUG
  PRECICE_DEBUG("Detected " << accessorGroups.size() << " groups");
  for (const AccessorGroup &group : accessorGroups) {
    PRECICE_DEBUG("Group " << group.id << ": name = " << group.name
                           << ", leaderRank = " << group.leaderRank
                           << ", size = " << group.size);
  }
#endif // NDEBUG
#endif // not PRECICE_NO_MPI
}

void Parallel::finalizeMPI()
{
  PRECICE_TRACE();
#ifndef PRECICE_NO_MPI
  if (_mpiInitializedByPrecice) {
    PRECICE_DEBUG("preCICE finalizes MPI");
    MPI_Finalize();
  }
#endif // not PRECICE_NO_MPI
  _isInitialized = false;
}

// void Parallel::clearGroups()
// {
//   _accessorGroups.clear();
//   _isSplit = false;
// }

int Parallel::getProcessRank()
{
  // Do not use TRACE or DEBUG here!
#ifndef PRECICE_NO_MPI
  if (!_isInitialized)
    return 0;

  return getGlobalCommState()->rank();
#else
  return 0;
#endif // not PRECICE_NO_MPI
}

int Parallel::getLocalProcessRank()
{
#ifndef PRECICE_NO_MPI
  if (!_isInitialized)
    return 0;

  return getLocalCommState()->rank();
#else
  return 0;
#endif // not PRECICE_NO_MPI
}

inline int Parallel::getCommunicatorSize()
{
  return getCommunicatorSize(_globalCommunicator);
}

int Parallel::getCommunicatorSize(Communicator comm)
{
  PRECICE_TRACE();
  int communicatorSize = 1;
#ifndef PRECICE_NO_MPI
  if (_isInitialized) {
    MPI_Comm_size(comm, &communicatorSize);
  }
#endif // not PRECICE_NO_MPI
  return communicatorSize;
}

// void Parallel::synchronizeProcesses()
// {
// #ifndef PRECICE_NO_MPI
//   PRECICE_TRACE();
//   PRECICE_ASSERT(_isInitialized);
//   MPI_Barrier(_globalCommunicator);
// #endif // not PRECICE_NO_MPI
// }
//
// void Parallel::synchronizeLocalProcesses()
// {
// #ifndef PRECICE_NO_MPI
//   PRECICE_TRACE();
//   PRECICE_ASSERT(_isInitialized && _isSplit);
//   MPI_Barrier(_localCommunicator);
// #endif // not PRECICE_NO_MPI
// }

// // @TODO The freeing and cleaning behaviour is weird and needs fixing
// void Parallel::setGlobalCommunicator(
//     Parallel::Communicator defaultCommunicator,
//     bool                   free)
// {
// #ifndef PRECICE_NO_MPI
//   PRECICE_TRACE();
//   if (free && _globalCommunicator != getCommunicatorWorld() && _globalCommunicator != MPI_COMM_SELF && _globalCommunicator != MPI_COMM_NULL) {
//     MPI_Comm_free(&_globalCommunicator);
//   }
//   _globalCommunicator = defaultCommunicator;
//   _localCommunicator  = _globalCommunicator;
//   if (free)
//     _accessorGroups.clear();
// #endif // not PRECICE_NO_MPI
// }

const Parallel::CommStatePtr Parallel::getGlobalCommState()
{
  PRECICE_TRACE();
  auto local = current();
  return local->parent ? local->parent->comm : MPI_COMM_WORLD;
}

const Parallel::CommStatePtr Parallel::getLocalCommState()
{
  PRECICE_TRACE();
  return current()->comm;
}

void Parallel::restrictCommunicator(int newSize)
{
  PRECICE_TRACE();
  PRECICE_ASSERT(newSize > 0, "Cannot restrict a communicator to nothing!");

#ifndef PRECICE_NO_MPI
  PRECICE_ASSERT(_isInitialized);
  auto       baseState = current();
  const auto size      = baseState->size();
  const auto rank      = baseState->rank();

  PRECICE_ASSERT(newSize > size, "Requested more ranks than the Communicator can provide!");

  // A single rank can use MPI_COMM_SELF, nothing else required
  if (newSize == 1) {
    // @todo verify if the barrier is required
    // PRECICE_DEBUG("Barrier");
    // MPI_Barrier(_globalCommunicator);
    if (rank == 0) {
      PRECICE_DEBUG("Restricted to COMM_SELF");
      pushState(CommState::self());
      return;
    } else {
      PRECICE_DEBUG("This rank remains unused after the restriction.");
      pushState(CommState::null());
      return;
    }
  }

  // If the requested size is the same as the capacity, then simply duplicate the comm
  if (newSize == size) {
    PRECICE_DEBUG("Restriction to capacity: duplicating Comm");
    MPI_Comm copiedComm;
    MPI_Comm_dup(baseState->comm, &copiesComm);
    pushState(CommState::fromComm(copiedComm));
    return;
  }

  // Create group, containing all processes of communicator
  MPI_Group currentGroup;
  MPI_Comm_group(baseState->comm, &currentGroup);

  // Prepare a ranks vector with 0,1,2,...,newSize-1
  std::vector<int> ranks(newSize);
  std::iota(ranks.begin(), ranks.end(), 0);

  // Create subgroup, containing processes contained in ranks
  PRECICE_DEBUG("Create restricted group");
  MPI_Group restrictedGroup;
  MPI_Group_incl(currentGroup, ranks.size(), ranks.data(), &restrictedGroup);

  // @todo check if it has to go back down to the other group free
  MPI_Group_free(&currentGroup);

  {
    int restrictedGroupSize = 0;
    MPI_Group_size(restrictedGroup, &restrictedGroupSize);
    PRECICE_ASSERT(restrictedGroupSize == newSize, "Group size differs: restriction failed!");
  }

  // Create communicator, containing process of restrictedGroup
  PRECICE_DEBUG("Create restricted comm using group");

  //  MPI_Comm restrictedCommunicator;
  MPI_Comm restrictedCommunicator = MPI_COMM_NULL;
  MPI_Comm_create(baseState->comm, restrictedGroup, &restrictedCommunicator);

  // @todo check if the we need a barrier here
  // PRECICE_DEBUG("Barrier");
  // MPI_Barrier(_globalCommunicator);
  PRECICE_DEBUG("Free restricted group");
  MPI_Group_free(&restrictedGroup);

  pushState(CommState::fromComm(restrictedCommunicator));
#endif // not PRECICE_NO_MPI
}

const std::vector<Parallel::AccessorGroup> &Parallel::getAccessorGroups()
{
  PRECICE_TRACE();
  PRECICE_ASSERT(_isInitialized);
  return _accessorGroups;
}
} // namespace utils
} // namespace precice

//#endif // not PRECICE_NO_MPI
