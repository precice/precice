#ifndef PRECICE_NO_MPI
#include "com/CommunicateBoundingBox.hpp"
#include "com/MPIDirectCommunication.hpp"
#include "mesh/Mesh.hpp"
#include "testing/Testing.hpp"
#include "testing/Fixtures.hpp"
#include "utils/Parallel.hpp"

using namespace precice;
using namespace precice::com;

BOOST_AUTO_TEST_SUITE(CommunicationTests)

BOOST_AUTO_TEST_SUITE(CommunicateBoundingBoxTests)

BOOST_FIXTURE_TEST_CASE(SendAndReceiveBoundingBox, testing::M2NFixture,
                       * testing::MinRanks(2)
                       * boost::unit_test::fixture<testing::MPICommRestrictFixture>(std::vector<int>({0, 1})))
{
  if (utils::Parallel::getCommunicatorSize() != 2)
    return;

  for (int dim = 2; dim <= 3; dim++) {
    mesh::Mesh::BoundingBox bb;

    for (int i=0; i < dim; i++) {
      bb.push_back(std::make_pair(i,i+1));
    }

    CommunicateBoundingBox comBB(m2n->getMasterCommunication());

    if (utils::Parallel::getProcessRank() == 0) {
      comBB.sendBoundingBox(bb, 0); // send to 0 as we communicate between both master ranks
    }
    else if (utils::Parallel::getProcessRank() == 1) {

      mesh::Mesh::BoundingBox bbCompare;
      for (int i=0; i < dim; i++) {
        bbCompare.push_back(std::make_pair(-1,-1));
      }

      comBB.receiveBoundingBox(bbCompare, 0);

      BOOST_TEST(bb==bbCompare);
    }
  }
}

//@todo: tests for all other methods

BOOST_AUTO_TEST_SUITE_END() // BB

BOOST_AUTO_TEST_SUITE_END() // Communication

#endif // not PRECICE_NO_MPI
