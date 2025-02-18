#ifndef PRECICE_NO_MPI

#include <vector>
#include "precice/impl/MeshContext.hpp"
#include "precice/impl/ParticipantImpl.hpp"
#include "precice/impl/ParticipantState.hpp"
#include "precice/precice.hpp"
#include "testing/Testing.hpp"

using namespace precice;

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(Whitebox)
/// Test reading of a full features coupling configuration file.
PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(TestConfigurationComsol)
{
  PRECICE_TEST();
  // Test configuration for accessor "Comsol"
  Participant interfaceComsol("Comsol", context.config(), 0, 1);
  BOOST_TEST(testing::WhiteboxAccessor::impl(interfaceComsol)._participants.size() == 2);

  impl::PtrParticipant comsol = testing::WhiteboxAccessor::impl(interfaceComsol)._participants.at(1);
  BOOST_TEST(comsol);
  BOOST_TEST(comsol->getName() == "Comsol");

  const auto &meshContexts = comsol->_meshContexts;
  BOOST_TEST(meshContexts.size() == 1);
  BOOST_TEST(meshContexts.count("PeanoNodes") == 0);
  BOOST_TEST(meshContexts.count("ComsolNodes") > 0);
  BOOST_TEST(meshContexts.at("ComsolNodes")->mesh->getName() == std::string("ComsolNodes"));
  BOOST_TEST(comsol->_usedMeshContexts.size() == 1);

  BOOST_TEST(interfaceComsol.getMeshDimensions("ComsolNodes") == 2);
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // Whitebox

#endif // PRECICE_NO_MPI
