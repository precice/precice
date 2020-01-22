#include "testing/Testing.hpp"
#include <cstdlib>
#include <exception>
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
  Par::setGlobalCommunicator(Par::getCommunicatorWorld());
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
  participants.emplace_back(std::move(participant));
}

void TestContext::setContextFrom(const Participant &p, int rank)
{
  this->name = p.name;
  this->size = p.size;
  this->rank = rank;
}

void TestContext::initialize(const Participants &participants)
{
  initializeMPI(participants);
  initializeEvents();
  initializePetsc();
}

void TestContext::initializeMPI(const TestContext::Participants &participants)
{
  const int globalRank = Par::getProcessRank();
  const int required   = std::accumulate(participants.begin(), participants.end(), 0, [](int total, const Participant &next) { return total + next.size; });
  if (required > Par::getCommunicatorSize()) {
    throw std::runtime_error{"This test requests more ranks than available"};
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

  // If there were multiple participant requested, we need to split the restricted comm
  if (participants.size() > 1) {
    int offset = 0;
    for (const auto &participant : participants) {
      const auto localRank = globalRank - offset;
      // Check if my global rank maps to this participant
      if (localRank < participant.size) {
        Par::splitCommunicator(participant.name);
        Par::setGlobalCommunicator(Par::getLocalCommunicator());
        setContextFrom(participant, localRank);
        return;
      }
      offset += participant.size;
    }
  }
}

void TestContext::initializePetsc()
{
  if (!invalid && _petsc) {
    precice::utils::Petsc::initialize(nullptr, nullptr);
  }
}

void TestContext::initializeEvents()
{
  if (!invalid && !_events) {
    precice::utils::EventRegistry::instance().initialize();
  }
}

} // namespace testing
} // namespace precice
