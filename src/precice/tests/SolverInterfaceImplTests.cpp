#include "testing/Testing.hpp"

#include "precice/impl/SolverInterfaceImpl.hpp"
#include "precice/SolverInterface.hpp"
#include "precice/config/Configuration.hpp"
#include "utils/MasterSlave.hpp"
#include "precice/impl/Participant.hpp"

#include <string>

using namespace precice;

struct SolverInterfaceTestFixture {

  std::string _pathToTests;

  void reset(){
     mesh::Mesh::resetGeometryIDsGlobally();
     mesh::Data::resetDataCount();
     impl::Participant::resetParticipantCount();
     utils::MasterSlave::reset();
  }

  SolverInterfaceTestFixture() {
    _pathToTests = testing::getPathToSources() + "/precice/tests/";
    reset();
  }
};

BOOST_AUTO_TEST_SUITE(PreciceTests)
BOOST_AUTO_TEST_SUITE(SolverInterfaceImpl)


BOOST_FIXTURE_TEST_CASE(MeshInfo, SolverInterfaceTestFixture)
{
  std::string filename = _pathToTests + "/configuration.xml";
  // Test configuration for accessor "Peano"
  SolverInterface interfacePeano ("Peano", 0, 1);
  config::Configuration config;
  xml::configure(config.getXMLTag(), filename);
  interfacePeano._impl->configure(config.getSolverInterfaceConfiguration());

  {
      auto id = interfacePeano.getMeshID("PeanoNodes");
      auto info = interfacePeano._impl->makeMeshInfo(id);
      std::ostringstream oss;
      oss << info;
      std::string expected = std::to_string(id) + ":\"PeanoNodes\"";
      BOOST_TEST(expected == oss.str());
  }
  {
      auto id = 1234; // Some non-existent mesh
      auto info = interfacePeano._impl->makeMeshInfo(id);
      std::ostringstream oss;
      oss << info;
      std::string expected = std::to_string(id) + ":<unknown>";
      BOOST_TEST(expected == oss.str());
  }
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
