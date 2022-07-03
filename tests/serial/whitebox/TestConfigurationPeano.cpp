#ifndef PRECICE_NO_MPI

#include <vector>
#include "precice/SolverInterface.hpp"
#include "precice/impl/MeshContext.hpp"
#include "precice/impl/Participant.hpp"
#include "precice/impl/SolverInterfaceImpl.hpp"
#include "testing/Testing.hpp"

using namespace precice;

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(Whitebox)
/// Test reading of a full features coupling configuration file.
BOOST_AUTO_TEST_CASE(TestConfigurationPeano)
{
  PRECICE_TEST(1_rank);
  // Test configuration for accessor "Peano"
  SolverInterface interfacePeano("Peano", context.config(), 0, 1);

  BOOST_TEST(testing::WhiteboxAccessor::impl(interfacePeano)._participants.size() == 2);
  BOOST_TEST(interfacePeano.getDimensions() == 2);

  impl::PtrParticipant peano = testing::WhiteboxAccessor::impl(interfacePeano)._participants.at(0);
  BOOST_TEST(peano);
  BOOST_TEST(peano->getName() == "Peano");

  std::vector<impl::MeshContext *> meshContexts = peano->_meshContexts;
  BOOST_TEST(meshContexts.size() == 2);
  BOOST_TEST(peano->_usedMeshContexts.size() == 2);

  BOOST_TEST(meshContexts.at(0)->mesh->getName() == std::string("PeanoNodes"));
  BOOST_TEST(meshContexts.at(1)->mesh->getName() == std::string("ComsolNodes"));
}

BOOST_AUTO_TEST_SUITE_END() // Whitebox
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // Integration

#endif // PRECICE_NO_MPI
