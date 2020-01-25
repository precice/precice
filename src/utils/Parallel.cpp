#include "Parallel.hpp"
#include <map>
#include "assertion.hpp"
#include "com/MPIDirectCommunication.hpp"

namespace precice {
namespace utils {

logging::Logger Parallel::_log("utils::Parallel");

Parallel::Communicator Parallel::_globalCommunicator = Parallel::getCommunicatorWorld();

Parallel::Communicator Parallel::_localCommunicator = MPI_COMM_NULL;

bool Parallel::_isInitialized           = false;
bool Parallel::_isSplit                 = false;
bool Parallel::_mpiInitializedByPrecice = false;

std::vector<Parallel::AccessorGroup> Parallel::_accessorGroups;

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
  if (_accessorGroups.empty()) {
    PRECICE_DEBUG("Exchange group information");
    //_accessorGroups.clear(); // Makes reinitialization possible
    std::map<std::string, int> groupMap; // map from names to group ID
    MPI_Comm                   globalComm = getGlobalCommunicator();
    int                        rank       = -1;
    MPI_Comm_rank(globalComm, &rank);
    int size = -1;
    MPI_Comm_size(globalComm, &size);

    bool severalGroups = false;

    if (size > 1) {
      MPI_Comm localComm = _isSplit ? getLocalCommunicator() : globalComm;
      PRECICE_INFO("<<<< Global " << getCommunicatorSize(globalComm) << "; Local " << getCommunicatorSize(localComm) << "; Same " << (localComm == globalComm));
      com::MPIDirectCommunication com{globalComm, localComm};
      if (rank == 0) {
        groupMap[groupName] = 0;
        AccessorGroup newGroup;
        newGroup.id         = 0;
        newGroup.size       = 1;
        newGroup.leaderRank = 0;
        newGroup.name       = groupName;
        _accessorGroups.push_back(newGroup);
        for (int i = 1; i < size; i++) {
          std::string name;
          com.receive(name, i); // Receive group name from all ranks
          if (groupMap.find(name) == groupMap.end()) {
            groupMap[name] = _accessorGroups.size();
            AccessorGroup newGroup;
            newGroup.id         = _accessorGroups.size();
            newGroup.size       = 1;
            newGroup.leaderRank = i;
            newGroup.name       = name;
            _accessorGroups.push_back(newGroup);
          } else {
            _accessorGroups[groupMap[name]].size++;
          }
        }
        auto groupCount = (int) _accessorGroups.size();
        MPI_Bcast(&groupCount, 1, MPI_INT, 0, globalComm);

        for (const AccessorGroup &group : _accessorGroups) {
          for (int i = 1; i < size; i++) {
            com.send(group.name, i);
            com.send(group.leaderRank, i);
            com.send(group.id, i);
            com.send(group.size, i);
          }
        }
        severalGroups = (_accessorGroups.size() > 1);
      } else { // rank != 0
        com.send(groupName, 0);
        int groupCount = -1;
        MPI_Bcast(&groupCount, 1, MPI_INT, 0, globalComm);
        severalGroups = groupCount > 1;
        for (int i = 0; i < groupCount; i++) {
          AccessorGroup newGroup;
          com.receive(newGroup.name, 0);
          com.receive(newGroup.leaderRank, 0);
          com.receive(newGroup.id, 0);
          com.receive(newGroup.size, 0);
          _accessorGroups.push_back(newGroup);
        }
      }

      if (severalGroups) {
        PRECICE_DEBUG("Split groups");
        int groupID = -1;
        for (auto group : _accessorGroups) {
          if (group.name == groupName) {
            groupID = group.id;
          }
        }
        PRECICE_ASSERT(groupID != -1);
        // Create a new communicator that contains only ranks of my group
        MPI_Comm_split(_globalCommunicator, groupID, getProcessRank(), &_localCommunicator);
      } else {
        _localCommunicator = _globalCommunicator;
      }

#ifndef NDEBUG
      PRECICE_DEBUG("Detected " << _accessorGroups.size() << " groups");
      for (const AccessorGroup &group : _accessorGroups) {
        PRECICE_DEBUG("Group " << group.id << ": name = " << group.name
                               << ", leaderRank = " << group.leaderRank
                               << ", size = " << group.size);
      }
#endif // NDEBUG
    }
  }
#endif // not PRECICE_NO_MPI

  _isSplit = true;
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

void Parallel::clearGroups()
{
  _accessorGroups.clear();
  _isSplit = false;
}

int Parallel::getProcessRank()
{
  // Do not use TRACE or DEBUG here!
  int processRank = 0;
#ifndef PRECICE_NO_MPI
  if (_isInitialized) {
    MPI_Comm_rank(_globalCommunicator, &processRank);
  }
#endif // not PRECICE_NO_MPI
  return processRank;
}

int Parallel::getLocalProcessRank()
{
  int processRank = 0;
#ifndef PRECICE_NO_MPI
  if (_isInitialized) {
    MPI_Comm_rank(_localCommunicator, &processRank);
  }
#endif
  return processRank;
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

void Parallel::synchronizeProcesses()
{
#ifndef PRECICE_NO_MPI
  PRECICE_TRACE();
  PRECICE_ASSERT(_isInitialized);
  MPI_Barrier(_globalCommunicator);
#endif // not PRECICE_NO_MPI
}

void Parallel::synchronizeLocalProcesses()
{
#ifndef PRECICE_NO_MPI
  PRECICE_TRACE();
  PRECICE_ASSERT(_isInitialized && _isSplit);
  MPI_Barrier(_localCommunicator);
#endif // not PRECICE_NO_MPI
}

// @TODO The freeing and cleaning behaviour is weird and needs fixing
void Parallel::setGlobalCommunicator(
    Parallel::Communicator defaultCommunicator,
    bool free)
{
#ifndef PRECICE_NO_MPI
  PRECICE_TRACE();
  if (free && _globalCommunicator != getCommunicatorWorld() && _globalCommunicator != MPI_COMM_SELF && _globalCommunicator != MPI_COMM_NULL) {
    MPI_Comm_free(&_globalCommunicator);
  }
  _globalCommunicator = defaultCommunicator;
  _localCommunicator  = _globalCommunicator;
  if (free) _accessorGroups.clear();
#endif // not PRECICE_NO_MPI
}

const Parallel::Communicator &Parallel::getGlobalCommunicator()
{
  PRECICE_TRACE();
  return _globalCommunicator;
}

const Parallel::Communicator &Parallel::getLocalCommunicator()
{
  PRECICE_TRACE();
  PRECICE_ASSERT(_isSplit);
  return _localCommunicator;
}

Parallel::Communicator Parallel::getRestrictedCommunicator(const std::vector<int> &ranks)
{
  PRECICE_TRACE();

#ifndef PRECICE_NO_MPI
  PRECICE_ASSERT(_isInitialized);
  PRECICE_ASSERT(not ranks.empty());
  PRECICE_ASSERT(ranks.size() <= static_cast<size_t>(getCommunicatorSize()));

  // Shortcut for running in isolation
  if (ranks.size() == 1) {
    PRECICE_DEBUG("Barrier");
    MPI_Barrier(_globalCommunicator);
    PRECICE_DEBUG("Restricted to COMM_SELF");
    return MPI_COMM_SELF;
  }

  // Create group, containing all processes of communicator
  MPI_Group currentGroup;
  MPI_Comm_group(_globalCommunicator, &currentGroup);
  std::vector<int> ranksArray(ranks.size());
#ifndef NDEBUG
  int communicatorSize = 0;
  MPI_Comm_size(_globalCommunicator, &communicatorSize);
  PRECICE_DEBUG("Comm size:" << communicatorSize);
  PRECICE_DEBUG("ranks size:" << (int) ranks.size());
//   PRECICE_ASSERT( (int)ranks.size() < communicatorSize );
#endif // Debug
  for (size_t i = 0; i < ranks.size(); i++) {
    PRECICE_DEBUG("Adding rank " << ranks[i]);
    PRECICE_ASSERT(ranks[i] >= 0);
    ranksArray[i] = ranks[i];
  }

  // Create subgroup, containing processes contained in ranks
  PRECICE_DEBUG("Restrict Group");
  MPI_Group restrictedGroup;
  MPI_Group_incl(currentGroup, ranksArray.size(), ranksArray.data(), &restrictedGroup);
#ifndef NDEBUG
  int restrictedGroupSize = 0;
  MPI_Group_size(restrictedGroup, &restrictedGroupSize);
  PRECICE_ASSERT(restrictedGroupSize > 0);
#endif

  // Create communicator, containing process of restrictedGroup
  PRECICE_DEBUG("Create Comm");

  //  MPI_Comm restrictedCommunicator;
  Communicator restrictedCommunicator = getCommunicatorWorld();
  MPI_Comm_create(_globalCommunicator, restrictedGroup, &restrictedCommunicator);
  PRECICE_DEBUG("Barrier");
  MPI_Barrier(_globalCommunicator);
  PRECICE_DEBUG("Free current group");
  MPI_Group_free(&currentGroup);
  PRECICE_DEBUG("Free restricted group");
  MPI_Group_free(&restrictedGroup);
  return restrictedCommunicator;
#else  // not PRECICE_NO_MPI
  return getCommunicatorWorld();
#endif // not PRECICE_NO_MPI
}

void Parallel::restrictGlobalCommunicator(const std::vector<int> &ranks)
{
  auto restrComm = getRestrictedCommunicator(ranks);
  if (std::find(ranks.begin(), ranks.end(), getProcessRank()) != ranks.end()) {
    setGlobalCommunicator(restrComm);
  } else {
    setGlobalCommunicator(MPI_COMM_NULL);
  }
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
