#include <testing/Testing.hpp>

#include <precice/config/Configuration.hpp>
#include <precice/config/ParticipantConfiguration.hpp>
#include <xml/ConfigParser.hpp>
#include "precice/impl/MeshContext.hpp"

using namespace precice;
using namespace precice::impl;

BOOST_AUTO_TEST_SUITE(PreciceTests)

BOOST_AUTO_TEST_SUITE(MeshContextTests)

using precice::impl::MeshContext;
using ParticipantName = std::string;
using MeshName        = std::string;
void meshContextTest(std::string configFile, std::map<ParticipantName, std::map<MeshName, precice::impl::MeshContext::Dynamicity>> meshSetup)
{
  BOOST_REQUIRE(!meshSetup.empty());
  config::Configuration     config;
  auto                      someParticipant = meshSetup.begin()->first;
  xml::ConfigurationContext cont{someParticipant, 0, 1};
  xml::configure(config.getXMLTag(),
                 cont,
                 configFile);
  auto participants = config.getParticipantConfiguration()->getParticipants();

  BOOST_REQUIRE(participants.size() == meshSetup.size());

  for (const auto &participant : participants) {
    auto pname = participant->getName();
    BOOST_TEST_CONTEXT(pname)
    {
      const auto &contextSetup = meshSetup.at(pname);

      const auto &contexts = participant->usedMeshContexts();
      for (const auto &context : contexts) {
        auto meshName = context->mesh->getName();
        BOOST_TEST_CONTEXT(meshName)
        {
          BOOST_TEST(context->dynamic == contextSetup.at(meshName));
        }
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(ExchangeStatic)
{
  PRECICE_TEST(1_rank);
  meshContextTest(
      context.prefix("meshcontext-static.xml"),
      {{"A", {{"MeshA", MeshContext::Dynamicity::No}}},
       {"B", {{"MeshA", MeshContext::Dynamicity::No}, {"MeshB", MeshContext::Dynamicity::No}}}});
}

BOOST_AUTO_TEST_CASE(ExchangeDynamic)
{
  PRECICE_TEST(1_rank);

  meshContextTest(
      context.prefix("meshcontext-dynamic.xml"),
      {{"A", {{"Dynamic", MeshContext::Dynamicity::Yes}}},
       {"B", {{"Dynamic", MeshContext::Dynamicity::Yes}, {"Transitive", MeshContext::Dynamicity::Transitively}}}});
}

BOOST_AUTO_TEST_CASE(ExchangeTransitive)
{
  PRECICE_TEST(1_rank);

  meshContextTest(
      context.prefix("meshcontext-transitive.xml"),
      {{"A", {{"Transitive", MeshContext::Dynamicity::Transitively}}},
       {"B", {{"Dynamic", MeshContext::Dynamicity::Yes}, {"Transitive", MeshContext::Dynamicity::Transitively}}}});
}

BOOST_AUTO_TEST_CASE(ExchangeDirectAccess)
{
  PRECICE_TEST(1_rank);

  meshContextTest(
      context.prefix("meshcontext-direct.xml"),
      {{"A", {{"Dynamic", MeshContext::Dynamicity::Yes}}},
       {"B", {{"Dynamic", MeshContext::Dynamicity::Yes}}}});
}

BOOST_AUTO_TEST_CASE(Partial)
{
  PRECICE_TEST(1_rank);

  meshContextTest(
      context.prefix("meshcontext-partial.xml"),
      {{"A", {{"A", MeshContext::Dynamicity::Transitively}, {"C", MeshContext::Dynamicity::No}}},
       {"B", {{"A", MeshContext::Dynamicity::Transitively}, {"B", MeshContext::Dynamicity::Yes}, {"C", MeshContext::Dynamicity::No}, {"D", MeshContext::Dynamicity::No}}}});
}

BOOST_AUTO_TEST_SUITE(Multi)

BOOST_AUTO_TEST_CASE(Static)
{
  PRECICE_TEST(1_rank);

  meshContextTest(
      context.prefix("meshcontext-multi-static.xml"),
      {
          {"SolverA", {
                          {"MeshA", MeshContext::Dynamicity::No},
                          {"MeshC", MeshContext::Dynamicity::No},
                      }},
          {"SolverB", {
                          {"MeshA", MeshContext::Dynamicity::No},
                          {"MeshC", MeshContext::Dynamicity::No},
                          {"MeshB1", MeshContext::Dynamicity::No},
                          {"MeshB2", MeshContext::Dynamicity::No},
                      }},
          {"SolverC", {
                          {"MeshC", MeshContext::Dynamicity::No},
                      }},
      });
}

BOOST_AUTO_TEST_CASE(Controller)
{
  PRECICE_TEST(1_rank);

  meshContextTest(
      context.prefix("meshcontext-multi-controller.xml"),
      {
          {"SolverA", {
                          {"MeshA", MeshContext::Dynamicity::Yes},
                          {"MeshC", MeshContext::Dynamicity::No},
                      }},
          {"SolverB", {
                          {"MeshA", MeshContext::Dynamicity::Yes},
                          {"MeshC", MeshContext::Dynamicity::No},
                          {"MeshB1", MeshContext::Dynamicity::Transitively},
                          {"MeshB2", MeshContext::Dynamicity::No},
                      }},
          {"SolverC", {
                          {"MeshC", MeshContext::Dynamicity::No},
                      }},
      });
}

BOOST_AUTO_TEST_CASE(Participant)
{
  PRECICE_TEST(1_rank);

  meshContextTest(
      context.prefix("meshcontext-multi-participant.xml"),
      {
          {"SolverA", {
                          {"MeshA", MeshContext::Dynamicity::No},
                          {"MeshC", MeshContext::Dynamicity::No},
                      }},
          {"SolverB", {
                          {"MeshA", MeshContext::Dynamicity::No},
                          {"MeshC", MeshContext::Dynamicity::Transitively},
                          {"MeshB1", MeshContext::Dynamicity::No},
                          {"MeshB2", MeshContext::Dynamicity::Yes},
                      }},
          {"SolverC", {
                          {"MeshC", MeshContext::Dynamicity::Transitively},
                      }},
      });
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
