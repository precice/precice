#include <algorithm>
#include <boost/filesystem.hpp>
#include <boost/filesystem/operations.hpp>
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
#include "logging/LogConfiguration.hpp"
#include "m2n/DistributedComFactory.hpp"
#include "m2n/GatherScatterComFactory.hpp"
#include "m2n/M2N.hpp"
#include "m2n/PointToPointComFactory.hpp"
#include "mesh/Data.hpp"
#include "precice/impl/Types.hpp"
#include "profiling/EventUtils.hpp"
#include "query/Index.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"
#include "utils/IntraComm.hpp"
#include "utils/Parallel.hpp"
#include "utils/Petsc.hpp"

namespace precice::testing {

using Par = utils::Parallel;

TestContext::~TestContext() noexcept
{
  if (!invalid && _petsc) {
    precice::utils::Petsc::finalize();
  }
  if (!invalid) {
    // Always clean up tests.
    precice::profiling::EventRegistry::instance().finalize();
  }
  if (!invalid && _initIntraComm) {
    utils::IntraComm::getCommunication() = nullptr;
    utils::IntraComm::reset();
  }

  // Reset communicators
  Par::resetCommState();
  Par::resetManagedMPI();
}

std::string TestContext::prefix(const std::string &filename) const
{
  boost::filesystem::path location{testing::getTestPath()};
  auto                    dir = location.parent_path();
  dir /= filename;
  return boost::filesystem::weakly_canonical(dir).string();
}

std::string TestContext::config() const
{
  auto testname = testing::getTestName();
  auto filepath = prefix(testname + ".xml");
  if (!boost::filesystem::is_regular_file(filepath)) {
    throw std::runtime_error("The requested config for \"" + testname + "\" does not exist at: " + filepath);
  }
  return filepath;
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

bool TestContext::isRank(Rank rank) const
{
  if (rank >= size) {
    throw std::runtime_error("The requested Rank does not exist!");
  }
  return this->rank == rank;
}

bool TestContext::isPrimary() const
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

void TestContext::handleOption(Participants &participants, ParticipantState participant)
{
  if (_simple) {
    std::terminate();
  }
  // @TODO add check if name already registered
  _names.push_back(participant.name);
  participants.emplace_back(std::move(participant));
}

void TestContext::setContextFrom(const ParticipantState &p)
{
  this->name           = p.name;
  this->size           = p.size;
  this->_initIntraComm = p.initIntraComm;
  this->_contextComm   = utils::Parallel::current();
  this->rank           = this->_contextComm->rank();
}

void TestContext::initialize(const Participants &participants)
{
  Par::Parallel::CommState::world()->synchronize();
  initializeMPI(participants);
  Par::Parallel::CommState::world()->synchronize();
  initializeIntraComm();
  initializeEvents();
  initializePetsc();
}

void TestContext::initializeMPI(const TestContext::Participants &participants)
{
  auto      baseComm   = Par::current();
  const int globalRank = baseComm->rank();
  const int available  = baseComm->size();

  // groups contain the accumulated sizes of previous groups
  std::vector<int> groups(participants.size());
  std::transform(participants.begin(), participants.end(), groups.begin(), [](const auto &p) { return p.size; });
  std::partial_sum(groups.begin(), groups.end(), groups.begin());

  // Check if there are enough ranks available
  auto required = groups.back();
  if (required > available) {
    throw std::runtime_error{"This test requests " + std::to_string(required) + " ranks, but there are only " + std::to_string(available) + " available"};
  }

  // Check if this rank isn't needed
  if (globalRank >= required) {
    Par::splitCommunicator(); // No group
    invalid = true;
    return;
  }

  // Find the participant this rank is assigned to
  auto position    = std::upper_bound(groups.begin(), groups.end(), globalRank);
  auto participant = std::distance(groups.begin(), position);
  Par::splitCommunicator(participant);
  setContextFrom(participants[participant]);
}

void TestContext::initializeIntraComm()
{
  if (invalid)
    return;

  // Establish a consistent state for all tests
  utils::IntraComm::configure(rank, size);
  utils::IntraComm::getCommunication().reset();
  logging::setMPIRank(rank);
  logging::setParticipant(name);

  if (!_initIntraComm || hasSize(1))
    return;

#ifndef PRECICE_NO_MPI
  precice::com::PtrCommunication intraComm = precice::com::PtrCommunication(new precice::com::MPIDirectCommunication());
#else
  precice::com::PtrCommunication intraComm = precice::com::PtrCommunication(new precice::com::SocketCommunication());
#endif

  intraComm->connectIntraComm(name, "", rank, size);

  utils::IntraComm::getCommunication() = std::move(intraComm);
}

void TestContext::initializeEvents()
{
  if (invalid) {
    return;
  }
  // Always initialize the events
  auto &er = precice::profiling::EventRegistry::instance();
  er.initialize(name, rank, size);
  if (_events) { // Enable them if they are requested
    er.setMode(precice::profiling::Mode::All);
    er.setDirectory("./precice-profiling");
  } else {
    er.setMode(precice::profiling::Mode::Off);
  }
  er.startBackend();
}

void TestContext::initializePetsc()
{
  if (!invalid && _petsc) {
    precice::utils::Petsc::initialize(_contextComm->comm);
  }
}

m2n::PtrM2N TestContext::connectPrimaryRanks(const std::string &acceptor, const std::string &requestor, const ConnectionOptions &options) const
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
  auto m2n = m2n::PtrM2N(new m2n::M2N(participantCom, distrFactory, options.useOnlyPrimaryCom, options.useTwoLevelInit));

  if (std::find(_names.begin(), _names.end(), acceptor) == _names.end()) {
    throw std::runtime_error{
        "Acceptor \"" + acceptor + "\" not defined in this context."};
  }
  if (std::find(_names.begin(), _names.end(), requestor) == _names.end()) {
    throw std::runtime_error{
        "Requestor \"" + requestor + "\" not defined in this context."};
  }

  if (isNamed(acceptor)) {
    m2n->acceptPrimaryRankConnection(acceptor, requestor);
  } else if (isNamed(requestor)) {
    m2n->requestPrimaryRankConnection(acceptor, requestor);
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

  if (_initIntraComm || _events || _petsc) {
    os << " Initialized: {";
    if (_initIntraComm)
      os << " IntraComm Communication ";
    if (_events)
      os << " Events";
    if (_petsc)
      os << " PETSc";
    os << '}';
  }
  return os.str();
}

} // namespace precice::testing
