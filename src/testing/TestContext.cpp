#include <algorithm>
#include <exception>
#include <filesystem>
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
#include "mapping/device/Ginkgo.hpp"
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

// TestSetup

void TestSetup::handleOption(testing::Require requirement)
{
  using testing::Require;
  switch (requirement) {
  case Require::PETSc:
    petsc  = true;
    events = true;
    break;
  case Require::Events:
    events = true;
    break;
  case Require::Ginkgo:
    ginkgo = true;
    events = true;
    break;
  default:
    std::terminate();
  }
}

void TestSetup::handleOption(ParticipantState participant)
{
  participants.emplace_back(participant);
}

void TestSetup::handleOption(Ranks ranks)
{
  participants.emplace_back("Unnamed"_on(ranks));
}

int TestSetup::totalRanks() const
{
  return std::accumulate(participants.begin(), participants.end(), 0, [](int i, const ParticipantState &ps) { return i + ps.size; });
}

// TestContext

TestContext::TestContext(TestSetup setup)
    : _setup(setup)
{
  for (const auto &p : setup.participants) {
    _names.emplace(p.name);
  }
  initialize(setup.participants);
}

TestContext::~TestContext() noexcept
{
  if (!invalid && _setup.petsc) {
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
  std::filesystem::path location{testing::getTestPath()};
  auto                  dir = location.parent_path();
  dir /= filename;
  return std::filesystem::weakly_canonical(dir).string();
}

std::string TestContext::config() const
{
  return prefix(testing::getTestName() + ".xml");
}

bool TestContext::hasSize(int size) const
{
  return this->size == size;
}

bool TestContext::isNamed(const std::string &name) const
{
  if (_names.count(name) == 0) {
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
  initializeGinkgo();
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
  if (_setup.events) { // Enable them if they are requested
    er.setMode(precice::profiling::Mode::All);
    er.setDirectory("./precice-profiling");
  } else {
    er.setMode(precice::profiling::Mode::Off);
  }
  er.startBackend();
}

void TestContext::initializePetsc()
{
  if (!invalid && _setup.petsc) {
    precice::utils::Petsc::initialize(_contextComm->comm);
  }
}

void TestContext::initializeGinkgo()
{
  if (!invalid && _setup.ginkgo) {
    int    argc = 0;
    char **argv;
#ifndef PRECICE_NO_GINKGO
    precice::device::Ginkgo::initialize(&argc, &argv);
#endif
  }
}

m2n::PtrM2N TestContext::connectPrimaryRanks(const std::string &acceptor, const std::string &requestor, const ConnectionOptions &options) const
{
  auto participantCom = com::PtrCommunication(new com::SocketCommunication());

  m2n::DistributedComFactory::SharedPointer distrFactory;
  switch (options.type) {
  case ConnectionType::GatherScatter:
    distrFactory = std::make_shared<m2n::GatherScatterComFactory>(participantCom);
    break;
  case ConnectionType::PointToPoint:
    distrFactory = std::make_shared<m2n::PointToPointComFactory>(com::PtrCommunicationFactory(new com::SocketCommunicationFactory()));
    break;
  default:
    throw std::runtime_error{"ConnectionType unknown"};
  };
  auto m2n = std::make_shared<m2n::M2N>(participantCom, distrFactory, options.useOnlyPrimaryCom, options.useTwoLevelInit);

  if (_names.count(acceptor) == 0) {
    throw std::runtime_error{
        "Acceptor \"" + acceptor + "\" not defined in this context."};
  }
  if (_names.count(requestor) == 0) {
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
  os << "Test context of " << testing::getFullTestName();
  if (name.empty()) {
    os << " is unnamed";
  } else {
    os << " represents \"" << name << '"';
  }
  os << " and runs on rank " << rank << " out of " << size << '.';

  if (_initIntraComm || _setup.events || _setup.petsc) {
    os << " Initialized: {";
    if (_initIntraComm)
      os << " IntraComm Communication ";
    if (_setup.events)
      os << " Events";
    if (_setup.petsc)
      os << " PETSc";
    os << '}';
  }
  return os.str();
}

} // namespace precice::testing
