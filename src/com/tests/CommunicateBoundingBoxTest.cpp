#ifndef PRECICE_NO_MPI
#include "com/CommunicateBoundingBox.hpp"
#include "mesh/Mesh.hpp"
#include "testing/Fixtures.hpp"
#include "testing/Testing.hpp"
#include "utils/Parallel.hpp"

using namespace precice;
using namespace precice::com;

BOOST_AUTO_TEST_SUITE(CommunicationTests)

BOOST_AUTO_TEST_SUITE(CommunicateBoundingBoxTests)

BOOST_FIXTURE_TEST_CASE(SendAndReceiveBoundingBox, testing::M2NFixture,
                        *testing::MinRanks(2) * boost::unit_test::fixture<testing::MPICommRestrictFixture>(std::vector<int>({0, 1})))
{
  if (utils::Parallel::getCommunicatorSize() != 2)
    return;

  for (int dim = 2; dim <= 3; dim++) {
    mesh::BoundingBox bb(dim);

    for (int i = 0; i < dim; i++) {
      bb.setBounds(i,i, i + 1);
    }

    CommunicateBoundingBox comBB(m2n->getMasterCommunication());

    if (utils::Parallel::getProcessRank() == 0) {
      comBB.sendBoundingBox(bb, 0);
    } else if (utils::Parallel::getProcessRank() == 1) {

      mesh::BoundingBox bbCompare(dim);
      for (int i = 0; i < dim; i++) {
        bbCompare.setBounds(i,-1, -1);
      }

      comBB.receiveBoundingBox(bbCompare, 0);

      BOOST_TEST(bb == bbCompare);
    }
  }
}

BOOST_FIXTURE_TEST_CASE(SendAndReceiveBoundingBoxMap, testing::M2NFixture,
                        *testing::MinRanks(2) * boost::unit_test::fixture<testing::MPICommRestrictFixture>(std::vector<int>({0, 1})))
{
  if (utils::Parallel::getCommunicatorSize() != 2)
    return;

  for (int dim = 2; dim <= 3; dim++) {
    mesh::BoundingBox    bb(dim);
    mesh::Mesh::BoundingBoxMap bbm;

    for (int rank = 0; rank < 3; rank++) {

      for (int i = 0; i < dim; i++) {
        bb.setBounds(i, rank * i, i + 1);
      }
      bbm[rank] = bb;
    }

    CommunicateBoundingBox comBB(m2n->getMasterCommunication());

    if (utils::Parallel::getProcessRank() == 0) {
      comBB.sendBoundingBoxMap(bbm, 0);
    } else if (utils::Parallel::getProcessRank() == 1) {

      mesh::BoundingBox    bbCompare(dim);
      mesh::Mesh::BoundingBoxMap bbmCompare;

      for (int rank = 0; rank < 3; rank++) {
        for (int i = 0; i < dim; i++) {
          bbCompare.setBounds(i,-1, -1);
        }
        bbmCompare[rank] = bbCompare;
      }

      comBB.receiveBoundingBoxMap(bbmCompare, 0);

      for (int rank = 0; rank < 3; rank++) {

        BOOST_TEST(bbm[rank] == bbmCompare[rank]);
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(BroadcastSendAndReceiveBoundingBoxMap,
                     *testing::OnSize(4) * boost::unit_test::fixture<testing::MasterComFixture>())
{

  // Build BB/BBMap to communicate

  mesh::BoundingBox    bb(3);
  mesh::Mesh::BoundingBoxMap bbm;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      bb.setBounds(j, i * j, i * (j + 1));
    }
    bbm[i] = bb;
  }

  CommunicateBoundingBox comBB(utils::MasterSlave::_communication);

  if (utils::Parallel::getProcessRank() == 0) {
    comBB.broadcastSendBoundingBoxMap(bbm);
  } else {

    mesh::BoundingBox    bbCompare(3);
    mesh::Mesh::BoundingBoxMap bbmCompare;

    for (int i = 0; i < 3; i++) {
      bbCompare.setBounds(i, -1, -1);
      bbmCompare[i] = bbCompare;
    }

    comBB.broadcastReceiveBoundingBoxMap(bbmCompare);

    for (int rank = 0; rank < 3; rank++) {

      BOOST_TEST(bbm[rank] == bbmCompare[rank]);
    }
  }
}

BOOST_FIXTURE_TEST_CASE(SendAndReceiveConnectionMap, testing::M2NFixture,
                        *testing::MinRanks(2) * boost::unit_test::fixture<testing::MPICommRestrictFixture>(std::vector<int>({0, 1})))
{
  if (utils::Parallel::getCommunicatorSize() != 2)
    return;

  std::vector<int>                fb;
  std::map<int, std::vector<int>> fbm;

  for (int rank = 0; rank < 3; rank++) {

    for (int i = 0; i < 3; i++) {
      fb.push_back(i + 1);
    }

    fbm[rank] = fb;
    fb.clear();
  }

  CommunicateBoundingBox comBB(m2n->getMasterCommunication());

  if (utils::Parallel::getProcessRank() == 0) {
    comBB.sendConnectionMap(fbm, 0);
  } else if (utils::Parallel::getProcessRank() == 1) {

    std::vector<int>                fbCompare;
    std::map<int, std::vector<int>> fbmCompare;

    for (int rank = 0; rank < 3; rank++) {

      for (int i = 0; i < 3; i++) {
        fbCompare.push_back(-1);
      }

      fbmCompare[rank] = fbCompare;
      fbCompare.clear();
    }

    comBB.receiveConnectionMap(fbmCompare, 0);

    for (int rank = 0; rank < 3; rank++) {

      BOOST_TEST(fbm[rank] == fbmCompare[rank]);
    }
  }
}

BOOST_AUTO_TEST_CASE(BroadcastSendAndReceiveConnectionMap,
                     *testing::OnSize(4) * boost::unit_test::fixture<testing::MasterComFixture>())
{

  std::vector<int>                fb;
  std::map<int, std::vector<int>> fbm;

  for (int rank = 0; rank < 3; rank++) {

    for (int i = 0; i < 3; i++) {
      fb.push_back(i + 1);
    }

    fbm[rank] = fb;
    fb.clear();
  }

  CommunicateBoundingBox comBB(utils::MasterSlave::_communication);

  if (utils::Parallel::getProcessRank() == 0) {
    comBB.broadcastSendConnectionMap(fbm);
  } else {

    std::vector<int>                fbCompare;
    std::map<int, std::vector<int>> fbmCompare;

    for (int rank = 0; rank < 3; rank++) {
      for (int i = 0; i < 3; i++) {
        fbCompare.push_back(-1);
      }
      fbmCompare[rank] = fbCompare;
      fbCompare.clear();
    }

    comBB.broadcastReceiveConnectionMap(fbmCompare);

    for (int rank = 0; rank < 3; rank++) {

      BOOST_TEST(fbm[rank] == fbmCompare[rank]);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END() // BB

BOOST_AUTO_TEST_SUITE_END() // Communication

#endif // not PRECICE_NO_MPI
