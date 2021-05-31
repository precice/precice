#ifndef PRECICE_NO_MPI
#include <map>
#include <memory>
#include <vector>
#include "com/CommunicateBoundingBox.hpp"
#include "com/SharedPointer.hpp"
#include "m2n/M2N.hpp"
#include "mesh/BoundingBox.hpp"
#include "mesh/Mesh.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"
#include "utils/MasterSlave.hpp"

using namespace precice;
using namespace precice::com;

BOOST_AUTO_TEST_SUITE(CommunicationTests)

BOOST_AUTO_TEST_SUITE(CommunicateBoundingBoxTests)

BOOST_AUTO_TEST_CASE(SendAndReceiveBoundingBox)
{
  PRECICE_TEST("A"_on(1_rank), "B"_on(1_rank), Require::Events);
  auto m2n = context.connectMasters("A", "B");

  for (int dim = 2; dim <= 3; dim++) {
    std::vector<double> bounds;
    for (int i = 0; i < dim; i++) {
      bounds.push_back(i);
      bounds.push_back(i + 1);
    }
    mesh::BoundingBox      bb(bounds);
    CommunicateBoundingBox comBB(m2n->getMasterCommunication());

    if (context.isNamed("A")) {
      comBB.sendBoundingBox(bb, 0);
    } else {
      BOOST_TEST(context.isNamed("B"));
      mesh::BoundingBox bbCompare(dim);

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
    mesh::Mesh::BoundingBoxMap bbm;

    for (int rank = 0; rank < 3; rank++) {
      std::vector<double> bounds;
      for (int i = 0; i < dim; i++) {
        bounds.push_back(rank * i);
        bounds.push_back(i + 1);
      }
      bbm.emplace(rank, mesh::BoundingBox(bounds));
    }

    CommunicateBoundingBox comBB(m2n->getMasterCommunication());

    if (context.isNamed("A")) {
      comBB.sendBoundingBoxMap(bbm, 0);
    } else {
      BOOST_TEST(context.isNamed("B"));

      mesh::BoundingBox          bbCompare(dim);
      mesh::Mesh::BoundingBoxMap bbmCompare;

      for (int i = 0; i < 3; i++) {
        bbmCompare.emplace(i, bbCompare);
      }

      comBB.receiveBoundingBoxMap(bbmCompare, 0);

      for (int rank = 0; rank < 3; rank++) {
        BOOST_TEST(bbm.at(rank) == bbmCompare.at(rank));
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(BroadcastSendAndReceiveBoundingBoxMap)
{
  PRECICE_TEST(""_on(4_ranks).setupMasterSlaves(), Require::Events);

  // Build BB/BBMap to communicate
  int                        dimension = 3;
  mesh::Mesh::BoundingBoxMap bbm;
  for (int rank = 0; rank < 3; rank++) {
    std::vector<double> bounds;
    for (int i = 0; i < dimension; i++) {
      bounds.push_back(rank * i);
      bounds.push_back(i + 1);
    }
    bbm.emplace(rank, mesh::BoundingBox(bounds));
  }

  CommunicateBoundingBox comBB(utils::MasterSlave::_communication);

  if (context.isMaster()) {
    comBB.broadcastSendBoundingBoxMap(bbm);
  } else {

    mesh::BoundingBox          bbCompare{dimension};
    mesh::Mesh::BoundingBoxMap bbmCompare;

    for (int i = 0; i < 3; i++) {
      bbmCompare.emplace(i, bbCompare);
    }
    comBB.broadcastReceiveBoundingBoxMap(bbmCompare);
    BOOST_TEST((int) bbmCompare.size() == 3);
    for (int rank = 0; rank < 3; rank++) {
      BOOST_TEST(bbm.at(rank) == bbmCompare.at(rank));
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

      BOOST_TEST(fbm.at(rank) == fbmCompare.at(rank));
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

      BOOST_TEST(fbm.at(rank) == fbmCompare.at(rank));
    }
  }
}

BOOST_AUTO_TEST_SUITE_END() // BB

BOOST_AUTO_TEST_SUITE_END() // Communication

#endif // not PRECICE_NO_MPI
