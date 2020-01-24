#include "testing/Testing.hpp"
#include <cstdlib>
#include <exception>
#include <string>
#include "com/MPIDirectCommunication.hpp"
#include "logging/LogMacros.hpp"
#include "utils/EventUtils.hpp"
#include "utils/Parallel.hpp"
#include "utils/Petsc.hpp"

namespace precice {
namespace testing {

std::string getPathToSources()
{
  precice::logging::Logger _log("testing");
  char *                   preciceRoot = std::getenv("PRECICE_ROOT");
  PRECICE_CHECK(preciceRoot != nullptr,
                "Environment variable PRECICE_ROOT has not been set. Please set it to the precice directory.");
  return std::string(preciceRoot) + "/src";
}

TestContext::~TestContext() noexcept
{
  if (!invalid && _petsc) {
    precice::utils::Petsc::finalize();
  }
  if (!invalid && _events) {
    precice::utils::EventRegistry::instance().finalize();
  }
  if (!invalid && _initMS) {
      utils::MasterSlave::_communication = nullptr;
      utils::MasterSlave::reset();
  }
  // Reset communicators
  Par::setGlobalCommunicator(Par::getCommunicatorWorld());
}

bool TestContext::hasSize(int size) const
{
  return this->size == size;
}

bool TestContext::isNamed(const std::string &name) const
{
  if (std::find(_names.begin(), _names.end(), name) == _names.end()) {
    throw std::runtime_error("The requested name \"" + name + "\" does not exist!");
  }
  return this->name == name;
}

bool TestContext::isRank(int rank) const
{
  if (rank >= size) {
    throw std::runtime_error("The requested Rank does not exist!");
  }
  return this->rank == rank;
}

bool TestContext::isMaster() const
{
  return isRank(0);
}

void TestContext::handleOption(Participants &, testing::Require requirement)
{
  using testing::Require;
  switch (requirement) {
  case Require::PETSc:
    _petsc  = true;
    _events = true;
    break;
  case Require::Events:
    _events = true;
    break;
  default:
    std::terminate();
  }
}

void TestContext::handleOption(Participants &participants, Participant participant)
{
  if (_simple) {
    std::terminate();
  }
  // @TODO add check if name already registered
  _names.push_back(participant.name);
  participants.emplace_back(std::move(participant));
}

void TestContext::setContextFrom(const Participant &p, int rank)
{
  this->name    = p.name;
  this->size    = p.size;
  this->rank    = rank;
  this->_initMS = p.initMS;
}

void TestContext::initialize(const Participants &participants)
{
  initializeMPI(participants);
  initializeMasterSlave();
  initializeEvents();
  initializePetsc();
}

void TestContext::initializeMPI(const TestContext::Participants &participants)
{
  const int globalRank = Par::getProcessRank();
  const int required   = std::accumulate(participants.begin(), participants.end(), 0, [](int total, const Participant &next) { return total + next.size; });
  const int available  = Par::getCommunicatorSize();
  if (required > available) {
    throw std::runtime_error{"This test requests " + std::to_string(required) +" ranks, but there are only "+std::to_string(available)+" available"};
  }

  // Restrict the communicator to the total required size
  Par::restrictGlobalCommunicator([required] {std::vector<int> v(required); std::iota(v.begin(), v.end(), 0); return v; }());

  // Mark all unnecessary ranks as invalid and return
  if (globalRank >= required) {
    invalid = true;
    Par::setGlobalCommunicator(MPI_COMM_NULL);
    return;
  }

  // If there was only a single participant requested, then update its info and we are done.
  if (participants.size() == 1) {
    auto &participant = participants.front();
    if (!invalid) {
      setContextFrom(participant, globalRank);
    }
    return;
  }

  // If there were multiple participants requested, we need to split the restricted comm
  if (participants.size() > 1) {
    int offset = 0;
    for (const auto &participant : participants) {
      const auto localRank = globalRank - offset;
      // Check if my global rank maps to this participant
      if (localRank < participant.size) {
        Par::splitCommunicator(participant.name);
        Par::setGlobalCommunicator(Par::getLocalCommunicator());
        Par::clearGroups();
        setContextFrom(participant, localRank);
        return;
      }
      offset += participant.size;
    }
  }
}

void TestContext::initializeMasterSlave()
{
  if (invalid || !_initMS) {
    return;
  }

  precice::com::PtrCommunication masterSlaveCom = precice::com::PtrCommunication(new precice::com::MPIDirectCommunication());
  utils::MasterSlave::_communication            = masterSlaveCom;
  utils::MasterSlave::configure(rank, size);

  const auto masterName = name + "Master";
  const auto slavesName = name + "Slaves";
  if (rank == 0) {
    masterSlaveCom->acceptConnection(masterName, slavesName, rank);
    masterSlaveCom->setRankOffset(1);
  } else {
    masterSlaveCom->requestConnection(masterName, slavesName, rank - 1, size - 1);
  }
}

void TestContext::initializeEvents()
{
  if (!invalid && _events) {
    precice::utils::EventRegistry::instance().initialize("precice-Tests", "", Par::getGlobalCommunicator());
  }
}

void TestContext::initializePetsc()
{
  if (!invalid && _petsc) {
    precice::utils::Petsc::initialize(nullptr, nullptr);
  }
}

} // namespace testing
} // namespace precice
