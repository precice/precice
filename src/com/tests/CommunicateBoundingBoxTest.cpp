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

BOOST_AUTO_TEST_SUITE(BoundingBoxTests)

BOOST_FIXTURE_TEST_CASE(TwoProcTestsWithM2NCommunication, testing::M2NFixture,
                       * testing::MinRanks(2)
                       * boost::unit_test::fixture<testing::MPICommRestrictFixture>(std::vector<int>({0, 1})))
{
  utils::Parallel::synchronizeProcesses();
  assertion(utils::Parallel::getCommunicatorSize() > 1);
  mesh::PropertyContainer::resetPropertyIDCounter();

  for (int dim = 2; dim <= 3; dim++) {
    // Build BB to communicate for rank0
    mesh::Mesh::BoundingBox BBToSend;

    if (utils::Parallel::getProcessRank() == 0) {
      for (int i=0; i < dim; i++) {
        BBToSend.push_back(std::make_pair(i,i+1));
      }      
    }

    if (utils::Parallel::getProcessRank() < 2) {

      com::PtrCommunication com(new com::MPIDirectCommunication());
      CommunicateBoundingBox comBB(com);
      
      if (utils::Parallel::getProcessRank() == 0) {

        comBB.sendBoundingBox(BBToSend, 0);        
      }
      
      else if (utils::Parallel::getProcessRank() == 1) {

        mesh::Mesh::BoundingBox BBToReceive, BBToCompare;
        for (int i=0; i < dim; i++) {
        BBToCompare.push_back(std::make_pair(i,i+1));
        BBToReceive.push_back(std::make_pair(0,0));
        }
      
        comBB.receiveBoundingBox(BBToReceive, 0, dim);
        BOOST_TEST(BBToReceive==BBToCompare);

        
      }
      com->closeConnection();
          
      utils::Parallel::setGlobalCommunicator(utils::Parallel::getCommunicatorWorld());
    }
  }
}
BOOST_AUTO_TEST_SUITE_END() // BB

BOOST_AUTO_TEST_SUITE_END() // Communication

#endif // not PRECICE_NO_MPI
