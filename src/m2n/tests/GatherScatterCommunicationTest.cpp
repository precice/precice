#ifndef PRECICE_NO_MPI

#include <Eigen/Core>
#include <memory>
#include <vector>
#include "com/SharedPointer.hpp"
#include "m2n/DistributedCommunication.hpp"
#include "m2n/M2N.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"

BOOST_AUTO_TEST_SUITE(M2NTests)

using namespace precice;
using namespace m2n;

BOOST_AUTO_TEST_CASE(GatherScatterTest)
{
  PRECICE_TEST("Part1"_on(1_rank), "Part2"_on(3_ranks).setupMasterSlaves(), Require::Events);
  auto m2n = context.connectMasters("Part1", "Part2");

  int             dimensions       = 2;
  int             numberOfVertices = 6;
  bool            flipNormals      = false;
  int             valueDimension   = 1;
  Eigen::VectorXd offset           = Eigen::VectorXd::Zero(dimensions);

  if (context.isNamed("Part1")) {
    mesh::PtrMesh pMesh(new mesh::Mesh("Mesh", dimensions, flipNormals, testing::nextMeshID()));
    m2n->createDistributedCommunication(pMesh);

    pMesh->setGlobalNumberOfVertices(numberOfVertices);
    pMesh->getVertexDistribution()[0].push_back(0);
    pMesh->getVertexDistribution()[0].push_back(1);
    pMesh->getVertexDistribution()[0].push_back(2);
    pMesh->getVertexDistribution()[0].push_back(3);
    pMesh->getVertexDistribution()[0].push_back(4);
    pMesh->getVertexDistribution()[0].push_back(5);

    m2n->acceptSlavesConnection("Part1", "Part2");
    Eigen::VectorXd values = Eigen::VectorXd::Zero(numberOfVertices);
    values << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0;
    m2n->send(values.data(), numberOfVertices, pMesh->getID(), valueDimension);
    m2n->receive(values.data(), numberOfVertices, pMesh->getID(),
                 valueDimension);
    BOOST_TEST(values(0) == 2.0);
    BOOST_TEST(values(1) == 4.0);
    BOOST_TEST(values(2) == 6.0);
    BOOST_TEST(values(3) == 16.0);
    BOOST_TEST(values(4) == 10.0);
    BOOST_TEST(values(5) == 12.0);

  } else {
    BOOST_TEST(context.isNamed("Part2"));
    mesh::PtrMesh pMesh(new mesh::Mesh("Mesh", dimensions, flipNormals, testing::nextMeshID()));
    m2n->createDistributedCommunication(pMesh);
    m2n->requestSlavesConnection("Part1", "Part2");

    if (context.isMaster()) {
      pMesh->setGlobalNumberOfVertices(numberOfVertices);
      pMesh->getVertexDistribution()[0].push_back(0);
      pMesh->getVertexDistribution()[0].push_back(1);
      pMesh->getVertexDistribution()[0].push_back(3);
      pMesh->getVertexDistribution()[2].push_back(2);
      pMesh->getVertexDistribution()[2].push_back(3);
      pMesh->getVertexDistribution()[2].push_back(4);
      pMesh->getVertexDistribution()[2].push_back(5);

      Eigen::Vector3d values(0.0, 0.0, 0.0);
      m2n->receive(values.data(), 3, pMesh->getID(), valueDimension);
      BOOST_TEST(values(0) == 1.0);
      BOOST_TEST(values(1) == 2.0);
      BOOST_TEST(values(2) == 4.0);
      values = values * 2;
      m2n->send(values.data(), 3, pMesh->getID(), valueDimension);
    } else if (context.isRank(1)) { // Slave1
      Eigen::VectorXd values;
      m2n->receive(values.data(), 0, pMesh->getID(), valueDimension);
      m2n->send(values.data(), 0, pMesh->getID(), valueDimension);
    } else {
      BOOST_TEST(context.isRank(2));
      Eigen::Vector4d values(0.0, 0.0, 0.0, 0.0);
      m2n->receive(values.data(), 4, pMesh->getID(), valueDimension);
      BOOST_TEST(values(0) == 3.0);
      BOOST_TEST(values(1) == 4.0);
      BOOST_TEST(values(2) == 5.0);
      BOOST_TEST(values(3) == 6.0);
      values = values * 2;
      m2n->send(values.data(), 4, pMesh->getID(), valueDimension);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()

#endif // PRECICE_NO_MPI
