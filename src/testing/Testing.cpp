#include "testing/Testing.hpp"
#include <cstdlib>
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
  if (!invalid) {
      precice::utils::Petsc::finalize();
      precice::utils::EventRegistry::instance().finalize();
  }
  Par::setGlobalCommunicator(Par::getCommunicatorWorld());
}

#if 0
void TestContext::handleOption(Participants &, EventsTag)

    void TestContext::handleOption(Participants &, EventsTag)
{
  PRECICE_ASSERT(!_simple);
  events = true;
}

void TestContext::handleOption(Participants &, PetscTag)
{
  PRECICE_ASSERT(!_simple);
  petsc = true;
}
#endif

void TestContext::handleOption(Participants &participants, Participant participant)
{
  PRECICE_ASSERT(!simple);
  participants.push_back(participant);
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
  const int required = std::accumulate(participants.begin(), participants.end(), 0, [](int total, const Participant &next) { return total + next.size; });
  if (required > Par::getCommunicatorSize()) {
    throw std::runtime_error{"This test requests more ranks than available"};
  }

  // Restrict the communicator to the total required size
  Par::restrictGlobalCommunicator([required] {std::vector<int> v(required); std::iota(v.begin(), v.end(), 0); return v; }());
  Par::setGlobalCommunicator(Par::getLocalCommunicator());

  const int globalRank = Par::getProcessRank();

  // Mark all unnecessary ranks as invalid
  if (globalRank >= required) {
    invalid = true;
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
  if (invalid) { // && !petsc) {
    return;
  }
  precice::utils::Petsc::initialize(nullptr, nullptr);
}

void TestContext::initializeEvents()
{
  if (invalid) { // && !events) {
    return;
  }
  precice::utils::EventRegistry::instance().initialize();
}

} // namespace testing
} // namespace precice
