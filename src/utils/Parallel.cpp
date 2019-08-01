#include "Parallel.hpp"
#include <map>
#include "assertion.hpp"
#include "com/MPIDirectCommunication.hpp"

namespace precice
{
namespace utils
{

logging::Logger Parallel::_log("utils::Parallel");

Parallel::Communicator Parallel::_globalCommunicator = Parallel::getCommunicatorWorld();

Parallel::Communicator Parallel::_localCommunicator = MPI_COMM_NULL;

bool Parallel::_isInitialized = false;
bool Parallel::_isSplit = false;
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
    int *argc,
    char ***argv)
{
#ifndef PRECICE_NO_MPI
  P_TRACE();
  int isMPIInitialized;
  MPI_Initialized(&isMPIInitialized);
  if (not isMPIInitialized) {
    P_DEBUG("Initialize MPI");
    _mpiInitializedByPrecice = true;
    MPI_Init(argc, argv);
  }
  _isInitialized = true;
#endif // not PRECICE_NO_MPI
}

void Parallel::splitCommunicator(const std::string &groupName)
{
#ifndef PRECICE_NO_MPI
  P_TRACE(groupName);

  // Exchange group information
  if (_accessorGroups.empty()) {
    P_DEBUG("Exchange group information");
    //_accessorGroups.clear(); // Makes reinitialization possible
    std::map<std::string, int> groupMap; // map from names to group ID
    MPI_Comm globalComm = getGlobalCommunicator();
    int rank = -1;
    MPI_Comm_rank(globalComm, &rank);
    int size = -1;
    MPI_Comm_size(globalComm, &size);

    bool severalGroups = false;

    if (size > 1) {
      com::MPIDirectCommunication com;
      if (rank == 0) {
        groupMap[groupName] = 0;
        AccessorGroup newGroup;
        newGroup.id = 0;
        newGroup.size = 1;
        newGroup.leaderRank = 0;
        newGroup.name = groupName;
        _accessorGroups.push_back(newGroup);
        for (int i = 1; i < size; i++) {
          std::string name;
          com.receive(name, i); // Receive group name from all ranks
          if (groupMap.find(name) == groupMap.end()) {
            groupMap[name] = _accessorGroups.size();
            AccessorGroup newGroup;
            newGroup.id = _accessorGroups.size();
            newGroup.size = 1;
            newGroup.leaderRank = i;
            newGroup.name = name;
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
        P_DEBUG("Split groups");
        int groupID = -1;
        for (auto group : _accessorGroups) {
          if (group.name == groupName) {
            groupID = group.id;
          }
        }
        P_ASSERT(groupID != -1);
        // Create a new communicator that contains only ranks of my group
        MPI_Comm_split(_globalCommunicator, groupID, getProcessRank(), &_localCommunicator);
      } else {
        _localCommunicator = _globalCommunicator;
      }

#ifndef NDEBUG
      P_DEBUG("Detected " << _accessorGroups.size() << " groups");
      for (const AccessorGroup &group : _accessorGroups) {
        P_DEBUG("Group " << group.id << ": name = " << group.name
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
  P_TRACE();
#ifndef PRECICE_NO_MPI
  if (_mpiInitializedByPrecice) {
    P_DEBUG("preCICE finalizes MPI");
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

int Parallel::getCommunicatorSize()
{
  P_TRACE();
  int communicatorSize = 1;
#ifndef PRECICE_NO_MPI
  if (_isInitialized) {
    MPI_Comm_size(_globalCommunicator, &communicatorSize);
  }
#endif // not PRECICE_NO_MPI
  return communicatorSize;
}

void Parallel::synchronizeProcesses()
{
#ifndef PRECICE_NO_MPI
  P_TRACE();
  P_ASSERT(_isInitialized);
  MPI_Barrier(_globalCommunicator);
#endif // not PRECICE_NO_MPI
}

void Parallel::synchronizeLocalProcesses()
{
#ifndef PRECICE_NO_MPI
  P_TRACE();
  P_ASSERT(_isInitialized && _isSplit);
  MPI_Barrier(_localCommunicator);
#endif // not PRECICE_NO_MPI
}

void Parallel::setGlobalCommunicator(
    Parallel::Communicator defaultCommunicator)
{
#ifndef PRECICE_NO_MPI
  P_TRACE();
  if (_globalCommunicator != getCommunicatorWorld()) {
    MPI_Comm_free(&_globalCommunicator);
  }
  _globalCommunicator = defaultCommunicator;
  _localCommunicator = _globalCommunicator;
  _accessorGroups.clear();
#endif // not PRECICE_NO_MPI
}

const Parallel::Communicator &Parallel::getGlobalCommunicator()
{
  P_TRACE();
  return _globalCommunicator;
}

const Parallel::Communicator &Parallel::getLocalCommunicator()
{
  P_TRACE();
  P_ASSERT(_isSplit);
  return _localCommunicator;
}

Parallel::Communicator Parallel::getRestrictedCommunicator(const std::vector<int> &ranks)
{
  P_TRACE();
  Communicator restrictedCommunicator = getCommunicatorWorld();
#ifndef PRECICE_NO_MPI
  P_ASSERT(_isInitialized);
  P_ASSERT(not ranks.empty());
  P_ASSERT(ranks.size() <= static_cast<size_t>(getCommunicatorSize()));
  // Create group, containing all processes of communicator
  MPI_Group currentGroup;
  MPI_Comm_group(_globalCommunicator, &currentGroup);
  std::vector<int> ranksArray(ranks.size());
#ifndef NDEBUG
  int communicatorSize = 0;
  MPI_Comm_size(_globalCommunicator, &communicatorSize);
  P_DEBUG("Comm size:" << communicatorSize);
  P_DEBUG("ranks size:" << (int) ranks.size());
//   P_ASSERT( (int)ranks.size() < communicatorSize );
#endif // Debug
  for (size_t i = 0; i < ranks.size(); i++) {
    P_DEBUG("Adding rank " << ranks[i]);
    P_ASSERT(ranks[i] >= 0);
    ranksArray[i] = ranks[i];
  }
  // Create subgroup, containing processes contained in ranks
  P_DEBUG("Restrict Group");
  MPI_Group restrictedGroup;
  MPI_Group_incl(currentGroup, ranks.size(), ranksArray.data(), &restrictedGroup);
#ifndef NDEBUG
  int restrictedGroupSize = 0;
  MPI_Group_size(restrictedGroup, &restrictedGroupSize);
  P_ASSERT(restrictedGroupSize > 0);
#endif
  // Create communicator, containing process of restrictedGroup
  P_DEBUG("Create Comm");
  //  MPI_Comm restrictedCommunicator;
  MPI_Comm_create(_globalCommunicator, restrictedGroup, &restrictedCommunicator);
  P_DEBUG("Barrier");
  MPI_Barrier(_globalCommunicator);
  P_DEBUG("Free current group");
  MPI_Group_free(&currentGroup);
  P_DEBUG("Free restricted group");
  MPI_Group_free(&restrictedGroup);
#endif // not PRECICE_NO_MPI
  return restrictedCommunicator;
}

void Parallel::restrictGlobalCommunicator(const std::vector<int> &ranks)
{
  auto restrComm = getRestrictedCommunicator(ranks);
  if (std::find(ranks.begin(), ranks.end(), getProcessRank()) != ranks.end())
    setGlobalCommunicator(restrComm);
}

const std::vector<Parallel::AccessorGroup> &Parallel::getAccessorGroups()
{
  P_TRACE();
  P_ASSERT(_isInitialized);
  return _accessorGroups;
}
}
} // precice, utils

//#endif // not PRECICE_NO_MPI
