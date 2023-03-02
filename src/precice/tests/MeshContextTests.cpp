#include <boost/test/tools/old/interface.hpp>
#include <testing/Testing.hpp>

#include <precice/config/Configuration.hpp>
#include <precice/config/ParticipantConfiguration.hpp>
#include <xml/ConfigParser.hpp>
#include "precice/impl/MeshContext.hpp"

using namespace precice;
using namespace precice::impl;

BOOST_AUTO_TEST_SUITE(PreciceTests)

BOOST_AUTO_TEST_SUITE(MeshContextTests)

BOOST_AUTO_TEST_CASE(StaticMesh)
{
  PRECICE_TEST(1_rank);

  config::Configuration     config;
  xml::ConfigurationContext cont{"A", 0, 1};
  xml::configure(config.getXMLTag(),
                 cont,
                 context.prefix("meshcontext-static.xml"));
  auto participants = config.getSolverInterfaceConfiguration().getParticipantConfiguration()->getParticipants();

  BOOST_REQUIRE(participants.size() == 2);

  for (const auto &participant : participants) {
    auto pname = participant->getName();
    BOOST_REQUIRE(pname == "A" || pname == "B");
    const auto &contexts = participant->usedMeshContexts();

    for (const auto &context : contexts) {
      auto meshName = context->mesh->getName();
      BOOST_REQUIRE(meshName == "MeshA" || meshName == "MeshB");
      using Dynamicity = precice::impl::MeshContext::Dynamicity;
      if (pname == "A") {
        if (meshName == "MeshA") {
          BOOST_TEST(context->provideMesh);
          BOOST_TEST((context->dynamic == Dynamicity::No));
        }
      }
      if (pname == "B") {
        if (meshName == "MeshA") {
          BOOST_TEST(!context->provideMesh);
          BOOST_TEST((context->dynamic == Dynamicity::No));
        }
        if (meshName == "MeshB") {
          BOOST_TEST(context->provideMesh);
          BOOST_TEST((context->dynamic == Dynamicity::No));
        }
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(DynamicMesh)
{
  PRECICE_TEST(1_rank);

  config::Configuration     config;
  xml::ConfigurationContext cont{"A", 0, 1};
  xml::configure(config.getXMLTag(),
                 cont,
                 context.prefix("meshcontext-dynamic.xml"));
  auto participants = config.getSolverInterfaceConfiguration().getParticipantConfiguration()->getParticipants();

  BOOST_REQUIRE(participants.size() == 2);

  for (const auto &participant : participants) {
    auto pname = participant->getName();
    BOOST_REQUIRE(pname == "A" || pname == "B");
    const auto &contexts = participant->usedMeshContexts();

    for (const auto &context : contexts) {
      auto meshName = context->mesh->getName();
      BOOST_REQUIRE(meshName == "Dynamic" || meshName == "Static");
      using Dynamicity = precice::impl::MeshContext::Dynamicity;
      if (pname == "A") {
        if (meshName == "Dynamic") {
          BOOST_TEST(context->provideMesh);
          BOOST_TEST((context->dynamic == Dynamicity::Yes));
        }
      }
      if (pname == "B") {
        if (meshName == "Dynamic") {
          BOOST_TEST(!context->provideMesh);
          BOOST_TEST((context->dynamic == Dynamicity::Yes));
        }
        if (meshName == "Static") {
          BOOST_TEST(context->provideMesh);
          BOOST_TEST((context->dynamic == Dynamicity::Transitively));
        }
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(Partial)
{
  PRECICE_TEST(1_rank);

  config::Configuration     config;
  xml::ConfigurationContext cont{"A", 0, 1};
  xml::configure(config.getXMLTag(),
                 cont,
                 context.prefix("meshcontext-partial.xml"));
  auto participants = config.getSolverInterfaceConfiguration().getParticipantConfiguration()->getParticipants();

  BOOST_REQUIRE(participants.size() == 2);

  for (const auto &participant : participants) {
    auto pname = participant->getName();
    BOOST_REQUIRE(pname == "A" || pname == "B");
    const auto &contexts = participant->usedMeshContexts();

    for (const auto &context : contexts) {
      auto meshName = context->mesh->getName();
      BOOST_REQUIRE(meshName == "MeshA" || meshName == "MeshB" || meshName == "MeshC" || meshName == "MeshD");
      using Dynamicity = precice::impl::MeshContext::Dynamicity;
      if (pname == "A") {
        if (meshName == "A") {
          BOOST_TEST(context->provideMesh);
          BOOST_TEST((context->dynamic == Dynamicity::Transitively));
        }
        if (meshName == "C") {
          BOOST_TEST(context->provideMesh);
          BOOST_TEST((context->dynamic == Dynamicity::No));
        }
      }
      if (pname == "B") {
        if (meshName == "A") {
          BOOST_TEST(!context->provideMesh);
          BOOST_TEST((context->dynamic == Dynamicity::Transitively));
        }
        if (meshName == "B") {
          BOOST_TEST(context->provideMesh);
          BOOST_TEST((context->dynamic == Dynamicity::Yes));
        }
        if (meshName == "C") {
          BOOST_TEST(!context->provideMesh);
          BOOST_TEST((context->dynamic == Dynamicity::No));
        }
        if (meshName == "D") {
          BOOST_TEST(context->provideMesh);
          BOOST_TEST((context->dynamic == Dynamicity::No));
        }
      }
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
