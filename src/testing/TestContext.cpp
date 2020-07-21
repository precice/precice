#include "testing/TestContext.hpp"
#include <algorithm>
#include <exception>
#include <memory>
#include <numeric>
#include <ostream>
#include <utility>
#include "com/Communication.hpp"
#include "com/MPIDirectCommunication.hpp"
#include "com/SharedPointer.hpp"
#include "com/SocketCommunication.hpp"
#include "com/SocketCommunicationFactory.hpp"
#include "m2n/DistributedComFactory.hpp"
#include "m2n/GatherScatterComFactory.hpp"
#include "m2n/M2N.hpp"
#include "m2n/PointToPointComFactory.hpp"
#include "mesh/Data.hpp"
#include "mesh/RTree.hpp"
#include "utils/EventUtils.hpp"
#include "utils/MasterSlave.hpp"
#include "utils/Parallel.hpp"
#include "utils/Petsc.hpp"

namespace precice {
namespace testing {

using Par = utils::Parallel;

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

  // Clear caches
  mesh::rtree::clear();

  // Reset static ids and counters
  mesh::Data::resetDataCount();

  // Reset communicators
  Par::resetCommState();
  Par::resetManagedMPI();
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
  this->name         = p.name;
  this->size         = p.size;
  this->rank         = rank;
  this->_initMS      = p.initMS;
  this->_contextComm = utils::Parallel::current();
}

void TestContext::initialize(const Participants &participants)
{
  Par::Parallel::CommState::world()->synchronize();
  initializeMPI(participants);
  Par::Parallel::CommState::world()->synchronize();
  initializeMasterSlave();
  initializeEvents();
  initializePetsc();
}

void TestContext::initializeMPI(const TestContext::Participants &participants)
{
  auto      baseComm   = Par::current();
  const int globalRank = baseComm->rank();
  const int available  = baseComm->size();
  const int required   = std::accumulate(participants.begin(), participants.end(), 0, [](int total, const Participant &next) { return total + next.size; });
  if (required > available) {
    throw std::runtime_error{"This test requests " + std::to_string(required) + " ranks, but there are only " + std::to_string(available) + " available"};
  }

  // Restrict the communicator to the total required size
  Par::restrictCommunicator(required);

  // Mark all unnecessary ranks as invalid and return
  if (globalRank >= required) {
    invalid = true;
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
        setContextFrom(participant, localRank);
        return;
      }
      offset += participant.size;
    }
  }
}

void TestContext::initializeMasterSlave()
{
  if (invalid)
    return;

  if (_initMS && hasSize(1)) {
    throw std::runtime_error{
        "Initializing a Master Slave communication does not make sense for serial participant \"" + name + "\"!"};
  }

  // Establish a consistent state for all tests
  utils::MasterSlave::configure(rank, size);

  if (!_initMS)
    return;

#ifndef PRECICE_NO_MPI
  precice::com::PtrCommunication masterSlaveCom = precice::com::PtrCommunication(new precice::com::MPIDirectCommunication());
#else
  precice::com::PtrCommunication masterSlaveCom = precice::com::PtrCommunication(new precice::com::SocketCommunication());
#endif

  masterSlaveCom->connectMasterSlaves(name, "", rank, size);

  utils::MasterSlave::_communication = std::move(masterSlaveCom);
}

void TestContext::initializeEvents()
{
  if (!invalid && _events) {
    precice::utils::EventRegistry::instance().initialize("precice-Tests", "", _contextComm->comm);
  }
}

void TestContext::initializePetsc()
{
  if (!invalid && _petsc) {
    precice::utils::Petsc::initialize(nullptr, nullptr, _contextComm->comm);
  }
}

m2n::PtrM2N TestContext::connectMasters(const std::string &acceptor, const std::string &requestor, const ConnectionOptions &options) const
{
  auto participantCom = com::PtrCommunication(new com::SocketCommunication());

  m2n::DistributedComFactory::SharedPointer distrFactory;
  switch (options.type) {
  case ConnectionType::GatherScatter:
    distrFactory.reset(new m2n::GatherScatterComFactory(participantCom));
    break;
  case ConnectionType::PointToPoint:
    distrFactory.reset(new m2n::PointToPointComFactory(com::PtrCommunicationFactory(new com::SocketCommunicationFactory())));
    break;
  default:
    throw std::runtime_error{"ConnectionType unknown"};
  };
  auto m2n = m2n::PtrM2N(new m2n::M2N(participantCom, distrFactory, options.useOnlyMasterCom, options.useTwoLevelInit));

  if (std::find(_names.begin(), _names.end(), acceptor) == _names.end()) {
    throw std::runtime_error{
        "Acceptor \"" + acceptor + "\" not defined in this context."};
  }
  if (std::find(_names.begin(), _names.end(), requestor) == _names.end()) {
    throw std::runtime_error{
        "Requestor \"" + requestor + "\" not defined in this context."};
  }

  if (isNamed(acceptor)) {
    m2n->acceptMasterConnection(acceptor, requestor);
  } else if (isNamed(requestor)) {
    m2n->requestMasterConnection(acceptor, requestor);
  } else {
    throw std::runtime_error{"You try to connect " + acceptor + " and " + requestor + ", but this context is named " + name};
  }
  return m2n;
}

std::string TestContext::describe() const
{
  if (invalid)
    return "This test context is invalid!";

  std::ostringstream os;
  os << "Test context";
  if (name.empty()) {
    os << " is unnamed";
  } else {
    os << " represents \"" << name << '"';
  }
  os << " and runs on rank " << rank << " out of " << size << '.';

  if (_initMS || _events || _petsc) {
    os << " Initialized: {";
    if (_initMS)
      os << " MasterSlave Communication ";
    if (_events)
      os << " Events";
    if (_petsc)
      os << " PETSc";
    os << '}';
  }
  return os.str();
}

} // namespace testing
} // namespace precice
