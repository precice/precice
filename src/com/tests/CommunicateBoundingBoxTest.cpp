#ifndef PRECICE_NO_MPI
#include "com/CommunicateBoundingBox.hpp"
#include "mesh/Mesh.hpp"
#include "testing/Testing.hpp"
#include "utils/Parallel.hpp"

using namespace precice;
using namespace precice::com;

BOOST_AUTO_TEST_SUITE(CommunicationTests)

BOOST_AUTO_TEST_SUITE(CommunicateBoundingBoxTests)

BOOST_AUTO_TEST_CASE(SendAndReceiveBoundingBox)
{
  PRECICE_TEST("A"_on(1_rank), "B"_on(1_rank), Require::Events);
  auto m2n = context.connectMasters("A", "B");

  for (int dim = 2; dim <= 3; dim++) {
    mesh::Mesh::BoundingBox bb;

    for (int i = 0; i < dim; i++) {
      bb.push_back(std::make_pair(i, i + 1));
    }

    CommunicateBoundingBox comBB(m2n->getMasterCommunication());

    if (context.isNamed("A")) {
      comBB.sendBoundingBox(bb, 0);
    } else {
      BOOST_TEST (context.isNamed("B"));
      mesh::Mesh::BoundingBox bbCompare;
      for (int i = 0; i < dim; i++) {
        bbCompare.push_back(std::make_pair(-1, -1));
      }

      comBB.receiveBoundingBox(bbCompare, 0);

      BOOST_TEST(bb == bbCompare);
    }
  }
}

BOOST_AUTO_TEST_CASE(SendAndReceiveBoundingBoxMap)
{
  PRECICE_TEST("A"_on(1_rank), "B"_on(1_rank), Require::Events);
  auto m2n = context.connectMasters("A", "B");

  for (int dim = 2; dim <= 3; dim++) {
    mesh::Mesh::BoundingBox    bb;
    mesh::Mesh::BoundingBoxMap bbm;

    for (int rank = 0; rank < 3; rank++) {

      for (int i = 0; i < dim; i++) {
        bb.push_back(std::make_pair(rank * i, i + 1));
      }

      bbm[rank] = bb;
      bb.clear();
    }

    CommunicateBoundingBox comBB(m2n->getMasterCommunication());

    if (context.isNamed("A")) {
      comBB.sendBoundingBoxMap(bbm, 0);
    } else {
      BOOST_TEST(context.isNamed("B"));

      mesh::Mesh::BoundingBox    bbCompare;
      mesh::Mesh::BoundingBoxMap bbmCompare;

      for (int rank = 0; rank < 3; rank++) {

        for (int i = 0; i < dim; i++) {
          bbCompare.push_back(std::make_pair(-1, -1));
        }

        bbmCompare[rank] = bbCompare;
        bbCompare.clear();
      }

      comBB.receiveBoundingBoxMap(bbmCompare, 0);

      for (int rank = 0; rank < 3; rank++) {

        BOOST_TEST(bbm[rank] == bbmCompare[rank]);
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(BroadcastSendAndReceiveBoundingBoxMap)
{
  PRECICE_TEST(""_on(4_ranks).setupMasterSlaves(), Require::Events);

  // Build BB/BBMap to communicate

  mesh::Mesh::BoundingBox    bb;
  mesh::Mesh::BoundingBoxMap bbm;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      bb.push_back(std::make_pair(i * j, i * (j + 1)));
    }
    bbm[i] = bb;
    bb.clear();
  }

  CommunicateBoundingBox comBB(utils::MasterSlave::_communication);

  if (context.isMaster()) {
    comBB.broadcastSendBoundingBoxMap(bbm);
  } else {

    mesh::Mesh::BoundingBox    bbCompare;
    mesh::Mesh::BoundingBoxMap bbmCompare;

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        bbCompare.push_back(std::make_pair(-1, -1));
      }
      bbmCompare[i] = bbCompare;
      bbCompare.clear();
    }

    comBB.broadcastReceiveBoundingBoxMap(bbmCompare);

    for (int rank = 0; rank < 3; rank++) {

      BOOST_TEST(bbm[rank] == bbmCompare[rank]);
    }
  }
}

BOOST_AUTO_TEST_CASE(SendAndReceiveConnectionMap)
{
  PRECICE_TEST("A"_on(1_rank), "B"_on(1_rank), Require::Events);
  auto m2n = context.connectMasters("A", "B");

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

  if (context.isNamed("A")) {
    comBB.sendConnectionMap(fbm, 0);
  } else if (context.isNamed("B")) {

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

BOOST_AUTO_TEST_CASE(BroadcastSendAndReceiveConnectionMap)
{
  PRECICE_TEST(""_on(4_ranks).setupMasterSlaves(), Require::Events);

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

  if (context.isMaster()) {
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
