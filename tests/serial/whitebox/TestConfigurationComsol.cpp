#ifndef PRECICE_NO_MPI

#include <vector>
#include "helpers.hpp"
#include "precice/SolverInterface.hpp"
#include "precice/impl/MeshContext.hpp"
#include "precice/impl/Participant.hpp"
#include "precice/impl/SolverInterfaceImpl.hpp"
#include "testing/Testing.hpp"

using namespace precice;

BOOST_AUTO_TEST_SUITE(PreciceTests)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(Whitebox)
BOOST_AUTO_TEST_CASE(TestConfigurationComsol)
{
  PRECICE_TEST(1_rank);
  // Test configuration for accessor "Comsol"
  SolverInterface interfaceComsol("Comsol", context.config(), 0, 1);
  BOOST_TEST(testing::WhiteboxAccessor::impl(interfaceComsol)._participants.size() == 2);
  BOOST_TEST(interfaceComsol.getDimensions() == 2);

  impl::PtrParticipant comsol = testing::WhiteboxAccessor::impl(interfaceComsol)._participants.at(1);
  BOOST_TEST(comsol);
  BOOST_TEST(comsol->getName() == "Comsol");

  std::vector<impl::MeshContext *> meshContexts = comsol->_meshContexts;
  BOOST_TEST(meshContexts.size() == 2);
  BOOST_TEST(meshContexts.at(0) == static_cast<void *>(nullptr));
  BOOST_TEST(meshContexts.at(1)->mesh->getName() == std::string("ComsolNodes"));
  BOOST_TEST(comsol->_usedMeshContexts.size() == 1);
}

BOOST_AUTO_TEST_SUITE_END() // PreciceTests
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // Whitebox

#endif // PRECICE_NO_MPI
