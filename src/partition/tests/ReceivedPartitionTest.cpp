#ifndef PRECICE_NO_MPI

#include <Eigen/Core>
#include <algorithm>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "com/Communication.hpp"
#include "com/Extra.hpp"
#include "com/SharedPointer.hpp"
#include "com/SocketCommunication.hpp"
#include "com/SocketCommunicationFactory.hpp"
#include "fixtures.hpp"
#include "m2n/DistributedComFactory.hpp"
#include "m2n/M2N.hpp"
#include "mapping/Mapping.hpp"
#include "mapping/NearestNeighborMapping.hpp"
#include "mapping/NearestProjectionMapping.hpp"
#include "mapping/PetRadialBasisFctMapping.hpp"
#include "mapping/SharedPointer.hpp"
#include "mapping/impl/BasisFunctions.hpp"
#include "math/constants.hpp"
#include "mesh/BoundingBox.hpp"
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Vertex.hpp"
#include "partition/Partition.hpp"
#include "partition/ProvidedPartition.hpp"
#include "partition/ReceivedPartition.hpp"
#include "precice/impl/Types.hpp"
#include "precice/impl/versions.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"
#include "utils/assertion.hpp"

namespace precice::mesh {
class Edge;
} // namespace precice::mesh

using namespace precice;
using namespace partition;
using precice::testing::TestContext;

BOOST_AUTO_TEST_SUITE(PartitionTests)
BOOST_AUTO_TEST_SUITE(ReceivedPartitionTests)

void createSolidzMesh2D(mesh::PtrMesh pSolidzMesh)
{
  int dimensions = 2;
  BOOST_TEST(pSolidzMesh);
  BOOST_TEST(pSolidzMesh->getDimensions() == dimensions);
  Eigen::VectorXd position(dimensions);

  position << 0.0, 0.0;
  mesh::Vertex &v1 = pSolidzMesh->createVertex(position);
  v1.setGlobalIndex(0);
  position << 0.0, 1.95;
  mesh::Vertex &v2 = pSolidzMesh->createVertex(position);
  v2.setGlobalIndex(1);
  position << 0.0, 2.1;
  mesh::Vertex &v3 = pSolidzMesh->createVertex(position);
  v3.setGlobalIndex(2);
  position << 0.0, 4.5;
  mesh::Vertex &v4 = pSolidzMesh->createVertex(position);
  v4.setGlobalIndex(3);
  position << 0.0, 5.95;
  mesh::Vertex &v5 = pSolidzMesh->createVertex(position);
  v5.setGlobalIndex(4);
  position << 0.0, 6.1;
  mesh::Vertex &v6 = pSolidzMesh->createVertex(position);
  v6.setGlobalIndex(5);
  pSolidzMesh->createEdge(v1, v2);
  pSolidzMesh->createEdge(v2, v3);
  pSolidzMesh->createEdge(v3, v4);
  pSolidzMesh->createEdge(v4, v5);
  pSolidzMesh->createEdge(v5, v6);
  pSolidzMesh->computeBoundingBox();
}

void createSolidzMesh2DSmall(mesh::PtrMesh pSolidzMesh)
{
  int dimensions = 2;
  BOOST_TEST(pSolidzMesh);
  BOOST_TEST(pSolidzMesh->getDimensions() == dimensions);
  Eigen::VectorXd position(dimensions);

  position << 0.0, 0.0;
  mesh::Vertex &v1 = pSolidzMesh->createVertex(position);
  position << 0.0, 3.0;
  mesh::Vertex &v2 = pSolidzMesh->createVertex(position);
  position << 0.0, 6.0;
  mesh::Vertex &v3 = pSolidzMesh->createVertex(position);
  pSolidzMesh->createEdge(v1, v2);
  pSolidzMesh->createEdge(v2, v3);
  pSolidzMesh->computeBoundingBox();
}

void createNastinMesh2D(mesh::PtrMesh pNastinMesh, Rank rank)
{
  int dimensions = 2;
  BOOST_TEST(pNastinMesh);
  BOOST_TEST(pNastinMesh->getDimensions() == dimensions);

  if (rank == 0) {
    Eigen::VectorXd position(dimensions);
    position << 0.0, 0.0;
    pNastinMesh->createVertex(position);
    position << 0.0, 2.0;
    pNastinMesh->createVertex(position);
  } else if (rank == 1) {
    // not at interface
  } else if (rank == 2) {
    Eigen::VectorXd position(dimensions);
    position << 0.0, 4.0;
    pNastinMesh->createVertex(position);
    position << 0.0, 6.0;
    pNastinMesh->createVertex(position);
  }
  pNastinMesh->computeBoundingBox();
}

void createNastinMesh2D2(mesh::PtrMesh pNastinMesh, Rank rank)
{
  int dimensions = 2;
  PRECICE_ASSERT(pNastinMesh.use_count() > 0);
  PRECICE_ASSERT(pNastinMesh->getDimensions() == dimensions);

  if (rank == 0) {
    Eigen::VectorXd position(dimensions);
    position << 0.10, 0.10;
    pNastinMesh->createVertex(position);
    position << 0.90, 0.90;
    pNastinMesh->createVertex(position);
  } else if (rank == 1) {
    // not at interface
  } else if (rank == 2) {

    Eigen::VectorXd position(dimensions);
    position << 2.1, 2.1;
    pNastinMesh->createVertex(position);
    position << 2.9, 2.9;
    pNastinMesh->createVertex(position);
  }
  pNastinMesh->computeBoundingBox();
}

void createSolidzMesh3D(mesh::PtrMesh pSolidzMesh)
{
  int             dimensions = 3;
  Eigen::VectorXd position(dimensions);
  BOOST_TEST(pSolidzMesh);
  BOOST_TEST(pSolidzMesh->getDimensions() == dimensions);

  position << 0.0, 0.0, -0.1;
  mesh::Vertex &v1 = pSolidzMesh->createVertex(position);
  v1.setGlobalIndex(0);
  position << -1.0, 0.0, 0.0;
  mesh::Vertex &v2 = pSolidzMesh->createVertex(position);
  v2.setGlobalIndex(1);
  position << 1.0, 0.0, 0.0;
  mesh::Vertex &v3 = pSolidzMesh->createVertex(position);
  v3.setGlobalIndex(2);
  position << 0.0, -1.0, 0.0;
  mesh::Vertex &v4 = pSolidzMesh->createVertex(position);
  v4.setGlobalIndex(3);
  position << 0.0, 1.0, 0.0;
  mesh::Vertex &v5 = pSolidzMesh->createVertex(position);
  v5.setGlobalIndex(4);
  mesh::Edge &e1 = pSolidzMesh->createEdge(v1, v2);
  mesh::Edge &e2 = pSolidzMesh->createEdge(v2, v4);
  mesh::Edge &e3 = pSolidzMesh->createEdge(v4, v1);
  mesh::Edge &e4 = pSolidzMesh->createEdge(v1, v3);
  mesh::Edge &e5 = pSolidzMesh->createEdge(v3, v5);
  mesh::Edge &e6 = pSolidzMesh->createEdge(v5, v1);
  pSolidzMesh->createTriangle(e1, e2, e3);
  pSolidzMesh->createTriangle(e4, e5, e6);
  pSolidzMesh->computeBoundingBox();
}

void createNastinMesh3D(mesh::PtrMesh pNastinMesh, Rank rank)
{
  int dimensions = 3;
  BOOST_TEST(pNastinMesh);
  BOOST_TEST(pNastinMesh->getDimensions() == dimensions);

  if (rank == 0) { //Primary
    Eigen::VectorXd position(dimensions);
    position << -1.0, -1.0, 0.0;
    pNastinMesh->createVertex(position);
    position << -0.75, -0.75, 0.5;
    pNastinMesh->createVertex(position);
  } else if (rank == 1) { //SecondaryRank1
    // secondary1 not at interface
  } else if (rank == 2) { //Secondary rank 2
    Eigen::VectorXd position(dimensions);
    position << 0.0, 0.0, -1.0;
    pNastinMesh->createVertex(position);
    position << 0.5, 0.5, 0.0;
    pNastinMesh->createVertex(position);
  }
  pNastinMesh->computeBoundingBox();
}

void createNastinMesh3D2(mesh::PtrMesh pNastinMesh, Rank rank)
{
  int dimensions = 3;
  PRECICE_ASSERT(pNastinMesh.use_count() > 0);
  PRECICE_ASSERT(pNastinMesh->getDimensions() == dimensions);

  if (rank == 0) {
    Eigen::VectorXd position(dimensions);
    position << 0.10, 0.10, 0.1;
    pNastinMesh->createVertex(position);
    position << 0.90, 0.90, 0.9;
    pNastinMesh->createVertex(position);
  } else if (rank == 1) {
    // not at interface
  } else if (rank == 2) {

    Eigen::VectorXd position(dimensions);
    position << 2.1, 2.1, 2.1;
    pNastinMesh->createVertex(position);
    position << 2.9, 2.9, 2.1;
    pNastinMesh->createVertex(position);
  }
  pNastinMesh->computeBoundingBox();
}

BOOST_AUTO_TEST_CASE(RePartitionNNBroadcastFilter2D)
{
  PRECICE_TEST("Solid"_on(1_rank), "Fluid"_on(3_ranks).setupIntraComm(), Require::Events);
  auto m2n = context.connectPrimaryRanks("Solid", "Fluid");

  int             dimensions = 2;
  Eigen::VectorXd offset     = Eigen::VectorXd::Zero(dimensions);

  if (context.isNamed("Solid")) { //SOLIDZ
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, testing::nextMeshID()));
    createSolidzMesh2D(pSolidzMesh);
    BOOST_TEST(pSolidzMesh->nVertices() == 6);
    ProvidedPartition part(pSolidzMesh);
    part.addM2N(m2n);
    part.communicate();
  } else {
    BOOST_TEST(context.isNamed("Fluid"));
    mesh::PtrMesh pNastinMesh(new mesh::Mesh("NastinMesh", dimensions, testing::nextMeshID()));
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, testing::nextMeshID()));

    mapping::PtrMapping boundingFromMapping = mapping::PtrMapping(
        new mapping::NearestNeighborMapping(mapping::Mapping::CONSISTENT, dimensions));
    mapping::PtrMapping boundingToMapping = mapping::PtrMapping(
        new mapping::NearestNeighborMapping(mapping::Mapping::CONSERVATIVE, dimensions));
    boundingFromMapping->setMeshes(pSolidzMesh, pNastinMesh);
    boundingToMapping->setMeshes(pNastinMesh, pSolidzMesh);

    createNastinMesh2D(pNastinMesh, context.rank);

    double safetyFactor = 0.1;

    ReceivedPartition part(pSolidzMesh, ReceivedPartition::ON_PRIMARY_RANK, safetyFactor);
    part.addM2N(m2n);
    part.addFromMapping(boundingFromMapping);
    part.addToMapping(boundingToMapping);
    part.communicate();
    part.compute();

    BOOST_TEST_CONTEXT(*pSolidzMesh)
    {
      // check if the sending and filtering worked right
      if (context.isPrimary()) { //Primary
        BOOST_TEST(pSolidzMesh->nVertices() == 2);
        BOOST_TEST(pSolidzMesh->edges().size() == 1);
      } else if (context.isRank(1)) { //SecondaryRank1
        BOOST_TEST(pSolidzMesh->nVertices() == 0);
        BOOST_TEST(pSolidzMesh->edges().size() == 0);
      } else if (context.isRank(2)) { //Secondary rank 2
        BOOST_TEST(pSolidzMesh->nVertices() == 2);
        BOOST_TEST(pSolidzMesh->edges().size() == 1);
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(RePartitionNNDoubleNode2D)
{
  PRECICE_TEST("Solid"_on(1_rank), "Fluid"_on(3_ranks).setupIntraComm(), Require::Events);
  auto m2n = context.connectPrimaryRanks("Solid", "Fluid");

  int             dimensions = 2;
  Eigen::VectorXd offset     = Eigen::VectorXd::Zero(dimensions);

  if (context.isNamed("Solid")) { //SOLIDZ
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, testing::nextMeshID()));
    createSolidzMesh2DSmall(pSolidzMesh);
    ProvidedPartition part(pSolidzMesh);
    part.addM2N(m2n);
    part.communicate();
  } else {
    BOOST_TEST(context.isNamed("Fluid"));
    mesh::PtrMesh pNastinMesh(new mesh::Mesh("NastinMesh", dimensions, testing::nextMeshID()));
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, testing::nextMeshID()));

    mapping::PtrMapping boundingFromMapping = mapping::PtrMapping(
        new mapping::NearestNeighborMapping(mapping::Mapping::CONSISTENT, dimensions));
    mapping::PtrMapping boundingToMapping = mapping::PtrMapping(
        new mapping::NearestNeighborMapping(mapping::Mapping::CONSERVATIVE, dimensions));
    boundingFromMapping->setMeshes(pSolidzMesh, pNastinMesh);
    boundingToMapping->setMeshes(pNastinMesh, pSolidzMesh);

    createNastinMesh2D(pNastinMesh, context.rank);

    double safetyFactor = 0.5;

    ReceivedPartition part(pSolidzMesh, ReceivedPartition::ON_SECONDARY_RANKS, safetyFactor);
    part.addM2N(m2n);
    part.addFromMapping(boundingFromMapping);
    part.addToMapping(boundingToMapping);
    part.communicate();
    part.compute();

    // check if the sending and filtering worked right
    if (context.isPrimary()) { //Primary
      BOOST_TEST(pSolidzMesh->nVertices() == 2);
      BOOST_TEST(pSolidzMesh->edges().size() == 1);
    } else if (context.isRank(1)) { //SecondaryRank1
      BOOST_TEST(pSolidzMesh->nVertices() == 0);
      BOOST_TEST(pSolidzMesh->edges().size() == 0);
    } else if (context.isRank(2)) { //Secondary rank 2
      BOOST_TEST(pSolidzMesh->nVertices() == 2);
      BOOST_TEST(pSolidzMesh->edges().size() == 1);
    }
  }
}

BOOST_AUTO_TEST_CASE(RePartitionNPPreFilterPostFilter2D)
{
  PRECICE_TEST("Solid"_on(1_rank), "Fluid"_on(3_ranks).setupIntraComm(), Require::Events);
  auto m2n = context.connectPrimaryRanks("Solid", "Fluid");

  int dimensions = 2;

  if (context.isNamed("Solid")) { //SOLIDZ
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, testing::nextMeshID()));
    createSolidzMesh2D(pSolidzMesh);
    ProvidedPartition part(pSolidzMesh);
    part.addM2N(m2n);
    part.communicate();
  } else {
    BOOST_TEST(context.isNamed("Fluid"));
    mesh::PtrMesh pNastinMesh(new mesh::Mesh("NastinMesh", dimensions, testing::nextMeshID()));
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, testing::nextMeshID()));

    mapping::PtrMapping boundingFromMapping = mapping::PtrMapping(
        new mapping::NearestProjectionMapping(mapping::Mapping::CONSISTENT, dimensions));
    mapping::PtrMapping boundingToMapping = mapping::PtrMapping(
        new mapping::NearestNeighborMapping(mapping::Mapping::CONSERVATIVE, dimensions));
    boundingFromMapping->setMeshes(pSolidzMesh, pNastinMesh);
    boundingToMapping->setMeshes(pNastinMesh, pSolidzMesh);

    createNastinMesh2D(pNastinMesh, context.rank);

    double            safetyFactor = 0.1;
    ReceivedPartition part(pSolidzMesh, ReceivedPartition::ON_PRIMARY_RANK, safetyFactor);
    part.addM2N(m2n);
    part.addFromMapping(boundingFromMapping);
    part.addToMapping(boundingToMapping);
    part.communicate();
    part.compute();

    BOOST_TEST_CONTEXT(*pSolidzMesh)
    {
      // check if the sending and filtering worked right
      if (context.isPrimary()) { //Primary
        BOOST_TEST(pSolidzMesh->nVertices() == 3);
        BOOST_TEST(pSolidzMesh->edges().size() == 2);
      } else if (context.isRank(1)) { //SecondaryRank1
        BOOST_TEST(pSolidzMesh->nVertices() == 0);
        BOOST_TEST(pSolidzMesh->edges().size() == 0);
      } else if (context.isRank(2)) { //Secondary rank 2
        BOOST_TEST(pSolidzMesh->nVertices() == 3);
        BOOST_TEST(pSolidzMesh->edges().size() == 2);
      }
    }
  }
}

#ifndef PRECICE_NO_PETSC
BOOST_AUTO_TEST_CASE(RePartitionRBFGlobal2D)
{
  PRECICE_TEST("Solid"_on(1_rank), "Fluid"_on(3_ranks).setupIntraComm(), Require::Events, Require::PETSc);
  auto m2n = context.connectPrimaryRanks("Solid", "Fluid");

  int dimensions = 2;

  if (context.isNamed("Solid")) { //SOLIDZ
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, testing::nextMeshID()));
    createSolidzMesh2D(pSolidzMesh);
    ProvidedPartition part(pSolidzMesh);
    part.addM2N(m2n);
    part.communicate();
  } else {
    BOOST_TEST(context.isNamed("Fluid"));
    mesh::PtrMesh pNastinMesh(new mesh::Mesh("NastinMesh", dimensions, testing::nextMeshID()));
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, testing::nextMeshID()));

    mapping::PtrMapping boundingFromMapping = mapping::PtrMapping(
        new mapping::PetRadialBasisFctMapping<mapping::ThinPlateSplines>(mapping::Mapping::CONSISTENT, dimensions,
                                                                         mapping::ThinPlateSplines(), {{false, false, false}}));
    mapping::PtrMapping boundingToMapping = mapping::PtrMapping(
        new mapping::NearestNeighborMapping(mapping::Mapping::CONSERVATIVE, dimensions));
    boundingFromMapping->setMeshes(pSolidzMesh, pNastinMesh);
    boundingToMapping->setMeshes(pNastinMesh, pSolidzMesh);

    createNastinMesh2D(pNastinMesh, context.rank);

    double            safetyFactor = 20.0;
    ReceivedPartition part(pSolidzMesh, ReceivedPartition::NO_FILTER, safetyFactor);
    part.addM2N(m2n);
    part.addFromMapping(boundingFromMapping);
    part.addToMapping(boundingToMapping);
    part.communicate();
    part.compute();

    BOOST_TEST_CONTEXT(*pSolidzMesh)
    {
      BOOST_TEST(pSolidzMesh->getVertexOffsets().size() == 3);
      BOOST_TEST(pSolidzMesh->getVertexOffsets().at(0) == 6);
      BOOST_TEST(pSolidzMesh->getVertexOffsets().at(1) == 6);
      BOOST_TEST(pSolidzMesh->getVertexOffsets().at(2) == 12);
      BOOST_TEST(pSolidzMesh->getGlobalNumberOfVertices() == 6);

      // check if the sending and filtering worked right
      if (context.isPrimary()) { //Primary
        BOOST_TEST(pSolidzMesh->nVertices() == 6);
        BOOST_TEST(pSolidzMesh->edges().size() == 5);
        BOOST_TEST(pSolidzMesh->vertex(0).isOwner() == true);
        BOOST_TEST(pSolidzMesh->vertex(1).isOwner() == true);
        BOOST_TEST(pSolidzMesh->vertex(2).isOwner() == true);
        BOOST_TEST(pSolidzMesh->vertex(3).isOwner() == false);
        BOOST_TEST(pSolidzMesh->vertex(4).isOwner() == false);
        BOOST_TEST(pSolidzMesh->vertex(5).isOwner() == false);
        BOOST_TEST(pSolidzMesh->vertex(0).getGlobalIndex() == 0);
        BOOST_TEST(pSolidzMesh->vertex(1).getGlobalIndex() == 1);
        BOOST_TEST(pSolidzMesh->vertex(2).getGlobalIndex() == 2);
        BOOST_TEST(pSolidzMesh->vertex(3).getGlobalIndex() == 3);
        BOOST_TEST(pSolidzMesh->vertex(4).getGlobalIndex() == 4);
        BOOST_TEST(pSolidzMesh->vertex(5).getGlobalIndex() == 5);
      } else if (context.isRank(1)) { //Secondary rank 2
        BOOST_TEST(pSolidzMesh->nVertices() == 0);
        BOOST_TEST(pSolidzMesh->edges().size() == 0);
      } else if (context.isRank(2)) { //Secondary rank 3
        BOOST_TEST(pSolidzMesh->nVertices() == 6);
        BOOST_TEST(pSolidzMesh->edges().size() == 5);
        BOOST_TEST(pSolidzMesh->vertex(0).isOwner() == false);
        BOOST_TEST(pSolidzMesh->vertex(1).isOwner() == false);
        BOOST_TEST(pSolidzMesh->vertex(2).isOwner() == false);
        BOOST_TEST(pSolidzMesh->vertex(3).isOwner() == true);
        BOOST_TEST(pSolidzMesh->vertex(4).isOwner() == true);
        BOOST_TEST(pSolidzMesh->vertex(5).isOwner() == true);
        BOOST_TEST(pSolidzMesh->vertex(0).getGlobalIndex() == 0);
        BOOST_TEST(pSolidzMesh->vertex(1).getGlobalIndex() == 1);
        BOOST_TEST(pSolidzMesh->vertex(2).getGlobalIndex() == 2);
        BOOST_TEST(pSolidzMesh->vertex(3).getGlobalIndex() == 3);
        BOOST_TEST(pSolidzMesh->vertex(4).getGlobalIndex() == 4);
        BOOST_TEST(pSolidzMesh->vertex(5).getGlobalIndex() == 5);
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(RePartitionRBFLocal2D1)
{
  PRECICE_TEST("Solid"_on(1_rank), "Fluid"_on(3_ranks).setupIntraComm(), Require::Events, Require::PETSc);
  auto m2n = context.connectPrimaryRanks("Solid", "Fluid");

  int dimensions = 2;

  if (context.isNamed("Solid")) { //SOLIDZ
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, testing::nextMeshID()));
    createSolidzMesh2D(pSolidzMesh);
    ProvidedPartition part(pSolidzMesh);
    part.addM2N(m2n);
    part.communicate();
  } else {
    BOOST_TEST(context.isNamed("Fluid"));
    mesh::PtrMesh pNastinMesh(new mesh::Mesh("NastinMesh", dimensions, testing::nextMeshID()));
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, testing::nextMeshID()));

    double supportRadius = 0.25;

    mapping::PtrMapping boundingFromMapping = mapping::PtrMapping(
        new mapping::PetRadialBasisFctMapping<mapping::CompactThinPlateSplinesC2>(mapping::Mapping::CONSISTENT, dimensions,
                                                                                  mapping::CompactThinPlateSplinesC2(supportRadius), {{false, false, false}}));
    mapping::PtrMapping boundingToMapping = mapping::PtrMapping(
        new mapping::NearestNeighborMapping(mapping::Mapping::CONSERVATIVE, dimensions));
    boundingFromMapping->setMeshes(pSolidzMesh, pNastinMesh);
    boundingToMapping->setMeshes(pNastinMesh, pSolidzMesh);

    createNastinMesh2D(pNastinMesh, context.rank);

    double            safetyFactor = 20.0;
    ReceivedPartition part(pSolidzMesh, ReceivedPartition::NO_FILTER, safetyFactor);
    part.addM2N(m2n);
    part.addFromMapping(boundingFromMapping);
    part.addToMapping(boundingToMapping);
    part.communicate();
    part.compute();

    BOOST_TEST_CONTEXT(*pSolidzMesh)
    {
      BOOST_TEST(pSolidzMesh->getVertexOffsets().size() == 3);
      BOOST_TEST(pSolidzMesh->getVertexOffsets().at(0) == 3);
      BOOST_TEST(pSolidzMesh->getVertexOffsets().at(1) == 3);
      BOOST_TEST(pSolidzMesh->getVertexOffsets().at(2) == 6);
      BOOST_TEST(pSolidzMesh->getGlobalNumberOfVertices() == 6);

      // check if the sending and filtering worked right
      if (context.isPrimary()) { //Primary
        BOOST_TEST(pSolidzMesh->nVertices() == 3);
        BOOST_TEST(pSolidzMesh->edges().size() == 2);
        BOOST_TEST(pSolidzMesh->vertex(0).isOwner() == true);
        BOOST_TEST(pSolidzMesh->vertex(1).isOwner() == true);
        BOOST_TEST(pSolidzMesh->vertex(2).isOwner() == true);
        BOOST_TEST(pSolidzMesh->vertex(0).getGlobalIndex() == 0);
        BOOST_TEST(pSolidzMesh->vertex(1).getGlobalIndex() == 1);
        BOOST_TEST(pSolidzMesh->vertex(2).getGlobalIndex() == 2);
      } else if (context.isRank(1)) { //Secondary rank 2
        BOOST_TEST(pSolidzMesh->nVertices() == 0);
        BOOST_TEST(pSolidzMesh->edges().size() == 0);
      } else if (context.isRank(2)) { //Secondary rank 3
        BOOST_TEST(pSolidzMesh->nVertices() == 3);
        BOOST_TEST(pSolidzMesh->edges().size() == 2);
        BOOST_TEST(pSolidzMesh->vertex(0).isOwner() == true);
        BOOST_TEST(pSolidzMesh->vertex(1).isOwner() == true);
        BOOST_TEST(pSolidzMesh->vertex(2).isOwner() == true);
        BOOST_TEST(pSolidzMesh->vertex(0).getGlobalIndex() == 3);
        BOOST_TEST(pSolidzMesh->vertex(1).getGlobalIndex() == 4);
        BOOST_TEST(pSolidzMesh->vertex(2).getGlobalIndex() == 5);
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(RePartitionRBFLocal2D2)
{
  PRECICE_TEST("Solid"_on(1_rank), "Fluid"_on(3_ranks).setupIntraComm(), Require::Events, Require::PETSc);
  auto m2n = context.connectPrimaryRanks("Solid", "Fluid");

  int dimensions = 2;

  if (context.isNamed("Solid")) { //SOLIDZ
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, testing::nextMeshID()));
    createSolidzMesh2D(pSolidzMesh);
    ProvidedPartition part(pSolidzMesh);
    part.addM2N(m2n);
    part.communicate();
  } else {
    BOOST_TEST(context.isNamed("Fluid"));
    mesh::PtrMesh pNastinMesh(new mesh::Mesh("NastinMesh", dimensions, testing::nextMeshID()));
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, testing::nextMeshID()));

    double supportRadius = 2.45;

    mapping::PtrMapping boundingFromMapping = mapping::PtrMapping(
        new mapping::PetRadialBasisFctMapping<mapping::CompactThinPlateSplinesC2>(mapping::Mapping::CONSISTENT, dimensions,
                                                                                  mapping::CompactThinPlateSplinesC2(supportRadius), {{false, false, false}}));
    mapping::PtrMapping boundingToMapping = mapping::PtrMapping(
        new mapping::NearestNeighborMapping(mapping::Mapping::CONSERVATIVE, dimensions));
    boundingFromMapping->setMeshes(pSolidzMesh, pNastinMesh);
    boundingToMapping->setMeshes(pNastinMesh, pSolidzMesh);

    createNastinMesh2D(pNastinMesh, context.rank);

    double            safetyFactor = 20.0;
    ReceivedPartition part(pSolidzMesh, ReceivedPartition::NO_FILTER, safetyFactor);
    part.addM2N(m2n);
    part.addFromMapping(boundingFromMapping);
    part.addToMapping(boundingToMapping);
    part.communicate();
    part.compute();

    BOOST_TEST_CONTEXT(*pSolidzMesh)
    {
      BOOST_TEST(pSolidzMesh->getVertexOffsets().size() == 3);
      BOOST_TEST(pSolidzMesh->getVertexOffsets().at(0) == 4);
      BOOST_TEST(pSolidzMesh->getVertexOffsets().at(1) == 4);
      BOOST_TEST(pSolidzMesh->getVertexOffsets().at(2) == 9);
      BOOST_TEST(pSolidzMesh->getGlobalNumberOfVertices() == 6);

      // check if the sending and filtering worked right
      if (context.isPrimary()) { //Primary
        BOOST_TEST(pSolidzMesh->nVertices() == 4);
        BOOST_TEST(pSolidzMesh->edges().size() == 3);
        BOOST_TEST(pSolidzMesh->vertex(0).isOwner() == true);
        BOOST_TEST(pSolidzMesh->vertex(1).isOwner() == true);
        BOOST_TEST(pSolidzMesh->vertex(2).isOwner() == true);
        BOOST_TEST(pSolidzMesh->vertex(3).isOwner() == false);
        BOOST_TEST(pSolidzMesh->vertex(0).getGlobalIndex() == 0);
        BOOST_TEST(pSolidzMesh->vertex(1).getGlobalIndex() == 1);
        BOOST_TEST(pSolidzMesh->vertex(2).getGlobalIndex() == 2);
        BOOST_TEST(pSolidzMesh->vertex(3).getGlobalIndex() == 3);
      } else if (context.isRank(1)) { //Secondary rank 2
        BOOST_TEST(pSolidzMesh->nVertices() == 0);
        BOOST_TEST(pSolidzMesh->edges().size() == 0);
      } else if (context.isRank(2)) { //Secondary rank 3
        BOOST_TEST(pSolidzMesh->nVertices() == 5);
        BOOST_TEST(pSolidzMesh->edges().size() == 4);
        BOOST_TEST(pSolidzMesh->vertex(0).isOwner() == false);
        BOOST_TEST(pSolidzMesh->vertex(1).isOwner() == false);
        BOOST_TEST(pSolidzMesh->vertex(2).isOwner() == true);
        BOOST_TEST(pSolidzMesh->vertex(3).isOwner() == true);
        BOOST_TEST(pSolidzMesh->vertex(4).isOwner() == true);
        BOOST_TEST(pSolidzMesh->vertex(0).getGlobalIndex() == 1);
        BOOST_TEST(pSolidzMesh->vertex(1).getGlobalIndex() == 2);
        BOOST_TEST(pSolidzMesh->vertex(2).getGlobalIndex() == 3);
        BOOST_TEST(pSolidzMesh->vertex(3).getGlobalIndex() == 4);
        BOOST_TEST(pSolidzMesh->vertex(4).getGlobalIndex() == 5);
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(RePartitionRBFLocal3D)
{
  PRECICE_TEST("Solid"_on(1_rank), "Fluid"_on(3_ranks).setupIntraComm(), Require::Events, Require::PETSc);
  auto m2n = context.connectPrimaryRanks("Solid", "Fluid");

  int dimensions = 3;

  if (context.isNamed("Solid")) { //SOLIDZ
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, testing::nextMeshID()));
    createSolidzMesh3D(pSolidzMesh);
    ProvidedPartition part(pSolidzMesh);
    part.addM2N(m2n);
    part.communicate();
  } else {
    BOOST_TEST(context.isNamed("Fluid"));
    mesh::PtrMesh pNastinMesh(new mesh::Mesh("NastinMesh", dimensions, testing::nextMeshID()));
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, testing::nextMeshID()));

    double supportRadius1 = 1.2;
    double supportRadius2 = 0.2;

    mapping::PtrMapping boundingFromMapping = mapping::PtrMapping(
        new mapping::PetRadialBasisFctMapping<mapping::CompactThinPlateSplinesC2>(mapping::Mapping::CONSISTENT, dimensions,
                                                                                  mapping::CompactThinPlateSplinesC2(supportRadius1), {{false, false, false}}));
    mapping::PtrMapping boundingToMapping = mapping::PtrMapping(
        new mapping::PetRadialBasisFctMapping<mapping::CompactThinPlateSplinesC2>(mapping::Mapping::CONSERVATIVE, dimensions,
                                                                                  mapping::CompactThinPlateSplinesC2(supportRadius2), {{false, false, false}}));
    boundingFromMapping->setMeshes(pSolidzMesh, pNastinMesh);
    boundingToMapping->setMeshes(pNastinMesh, pSolidzMesh);

    createNastinMesh3D(pNastinMesh, context.rank);

    double            safetyFactor = 20.0;
    ReceivedPartition part(pSolidzMesh, ReceivedPartition::NO_FILTER, safetyFactor);
    part.addM2N(m2n);
    part.addFromMapping(boundingFromMapping);
    part.addToMapping(boundingToMapping);
    part.communicate();
    part.compute();

    BOOST_TEST_CONTEXT(*pSolidzMesh)
    {
      BOOST_TEST(pSolidzMesh->getVertexOffsets().size() == 3);
      BOOST_TEST(pSolidzMesh->getVertexOffsets().at(0) == 5);
      BOOST_TEST(pSolidzMesh->getVertexOffsets().at(1) == 5);
      BOOST_TEST(pSolidzMesh->getVertexOffsets().at(2) == 10);
      BOOST_TEST(pSolidzMesh->getGlobalNumberOfVertices() == 5);

      // check if the sending and filtering worked right
      if (context.isPrimary()) { //Primary
        BOOST_TEST(pSolidzMesh->nVertices() == 5);
        BOOST_TEST(pSolidzMesh->edges().size() == 6);
        BOOST_TEST(pSolidzMesh->triangles().size() == 2);
        BOOST_TEST(pSolidzMesh->vertex(0).isOwner() == true);
        BOOST_TEST(pSolidzMesh->vertex(1).isOwner() == true);
        BOOST_TEST(pSolidzMesh->vertex(2).isOwner() == false);
        BOOST_TEST(pSolidzMesh->vertex(3).isOwner() == false);
        BOOST_TEST(pSolidzMesh->vertex(4).isOwner() == false);
        BOOST_TEST(pSolidzMesh->vertex(0).getGlobalIndex() == 0);
        BOOST_TEST(pSolidzMesh->vertex(1).getGlobalIndex() == 1);
        BOOST_TEST(pSolidzMesh->vertex(2).getGlobalIndex() == 2);
        BOOST_TEST(pSolidzMesh->vertex(3).getGlobalIndex() == 3);
        BOOST_TEST(pSolidzMesh->vertex(4).getGlobalIndex() == 4);
      } else if (context.isRank(1)) { //Secondary rank 2
        BOOST_TEST(pSolidzMesh->nVertices() == 0);
        BOOST_TEST(pSolidzMesh->edges().size() == 0);
        BOOST_TEST(pSolidzMesh->triangles().size() == 0);
      } else if (context.isRank(2)) { //Secondary rank 3
        BOOST_TEST(pSolidzMesh->nVertices() == 5);
        BOOST_TEST(pSolidzMesh->edges().size() == 6);
        BOOST_TEST(pSolidzMesh->triangles().size() == 2);
        BOOST_TEST(pSolidzMesh->vertex(0).isOwner() == false);
        BOOST_TEST(pSolidzMesh->vertex(1).isOwner() == false);
        BOOST_TEST(pSolidzMesh->vertex(2).isOwner() == true);
        BOOST_TEST(pSolidzMesh->vertex(3).isOwner() == true);
        BOOST_TEST(pSolidzMesh->vertex(4).isOwner() == true);
        BOOST_TEST(pSolidzMesh->vertex(0).getGlobalIndex() == 0);
        BOOST_TEST(pSolidzMesh->vertex(1).getGlobalIndex() == 1);
        BOOST_TEST(pSolidzMesh->vertex(2).getGlobalIndex() == 2);
        BOOST_TEST(pSolidzMesh->vertex(3).getGlobalIndex() == 3);
        BOOST_TEST(pSolidzMesh->vertex(4).getGlobalIndex() == 4);
      }
    }
  }
}

#endif // PRECICE_NO_PETSC

BOOST_AUTO_TEST_CASE(RePartitionNPBroadcastFilter3D)
{
  PRECICE_TEST("Fluid"_on(3_ranks).setupIntraComm(), "Solid"_on(1_rank), Require::Events);
  auto m2n = context.connectPrimaryRanks("Solid", "Fluid");

  int dimensions = 3;

  if (context.isNamed("Solid")) { //SOLIDZ
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, testing::nextMeshID()));
    createSolidzMesh3D(pSolidzMesh);
    ProvidedPartition part(pSolidzMesh);
    part.addM2N(m2n);
    part.communicate();
  } else {
    BOOST_TEST(context.isNamed("Fluid"));
    mesh::PtrMesh pNastinMesh(new mesh::Mesh("NastinMesh", dimensions, testing::nextMeshID()));
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, testing::nextMeshID()));

    mapping::PtrMapping boundingFromMapping = mapping::PtrMapping(
        new mapping::NearestProjectionMapping(mapping::Mapping::CONSISTENT, dimensions));
    mapping::PtrMapping boundingToMapping = mapping::PtrMapping(
        new mapping::NearestNeighborMapping(mapping::Mapping::CONSERVATIVE, dimensions));
    boundingFromMapping->setMeshes(pSolidzMesh, pNastinMesh);
    boundingToMapping->setMeshes(pNastinMesh, pSolidzMesh);

    createNastinMesh3D(pNastinMesh, context.rank);

    double            safetyFactor = 20.0;
    ReceivedPartition part(pSolidzMesh, ReceivedPartition::ON_PRIMARY_RANK, safetyFactor);
    part.addM2N(m2n);
    part.addFromMapping(boundingFromMapping);
    part.addToMapping(boundingToMapping);
    part.communicate();
    part.compute();

    // check if the sending and filtering worked right
    if (context.isPrimary()) { //Primary
      BOOST_TEST(pSolidzMesh->nVertices() == 2);
      BOOST_TEST(pSolidzMesh->edges().size() == 1);
      BOOST_TEST(pSolidzMesh->triangles().size() == 0);
    } else if (context.isRank(1)) { //SecondaryRank1
      BOOST_TEST(pSolidzMesh->nVertices() == 0);
      BOOST_TEST(pSolidzMesh->edges().size() == 0);
      BOOST_TEST(pSolidzMesh->triangles().size() == 0);
    } else if (context.isRank(2)) { //Secondary rank 2
      BOOST_TEST(pSolidzMesh->nVertices() == 3);
      BOOST_TEST(pSolidzMesh->edges().size() == 3);
      BOOST_TEST(pSolidzMesh->triangles().size() == 1);
    }
  }
}

BOOST_AUTO_TEST_CASE(TestRepartitionAndDistribution2D)
{
  PRECICE_TEST("Solid"_on(1_rank), "Fluid"_on(3_ranks).setupIntraComm(), Require::Events);
  auto m2n = context.connectPrimaryRanks("Solid", "Fluid");

  int dimensions = 2;

  if (context.isNamed("Solid")) { //SOLIDZ
    mesh::PtrMesh pMesh(new mesh::Mesh("SolidzMesh", dimensions, testing::nextMeshID()));

    Eigen::VectorXd position(dimensions);
    position << 0.0, 0.0;
    pMesh->createVertex(position);
    position << 1.0, 0.0;
    pMesh->createVertex(position);
    position << 2.0, 0.0;
    pMesh->createVertex(position);

    pMesh->computeBoundingBox();

    ProvidedPartition part(pMesh);
    part.addM2N(m2n);
    part.communicate();

  } else {
    BOOST_TEST(context.isNamed("Fluid"));
    mesh::PtrMesh pMesh(new mesh::Mesh("NastinMesh", dimensions, testing::nextMeshID()));
    mesh::PtrMesh pOtherMesh(new mesh::Mesh("SolidzMesh", dimensions, testing::nextMeshID()));

    mapping::PtrMapping boundingFromMapping = mapping::PtrMapping(
        new mapping::NearestNeighborMapping(mapping::Mapping::CONSISTENT, dimensions));
    boundingFromMapping->setMeshes(pMesh, pOtherMesh);

    if (context.isPrimary()) { //Primary
      Eigen::VectorXd position(dimensions);
      position << 0.0, 0.0;
      pOtherMesh->createVertex(position);
      position << 0.8, 0.0;
      pOtherMesh->createVertex(position);
    } else if (context.isRank(1)) { //Secondary rank 2
      Eigen::VectorXd position(dimensions);
      position << 1.0, 0.0;
      pOtherMesh->createVertex(position);
      position << 1.2, 0.0;
      pOtherMesh->createVertex(position);
    } else if (context.isRank(2)) { //Secondary rank 3
      // no vertices
    }

    pOtherMesh->computeBoundingBox();

    double            safetyFactor = 20.0;
    ReceivedPartition part(pMesh, ReceivedPartition::ON_PRIMARY_RANK, safetyFactor);
    part.addM2N(m2n);
    part.addFromMapping(boundingFromMapping);
    part.communicate();
    part.compute();

    BOOST_TEST(pMesh->getVertexOffsets().size() == 3);
    BOOST_TEST(pMesh->getVertexOffsets().at(0) == 2);
    BOOST_TEST(pMesh->getVertexOffsets().at(1) == 3);
    BOOST_TEST(pMesh->getVertexOffsets().at(2) == 3);

    if (context.isPrimary()) { //Primary
      BOOST_TEST(pMesh->getVertexDistribution().at(0).size() == 2);
      BOOST_TEST(pMesh->getVertexDistribution().at(1).size() == 1);
      BOOST_TEST(pMesh->getVertexDistribution().at(2).size() == 0);
      BOOST_TEST(pMesh->getVertexDistribution().at(0).at(0) == 0);
      BOOST_TEST(pMesh->getVertexDistribution().at(0).at(1) == 1);
      BOOST_TEST(pMesh->getVertexDistribution().at(1).at(0) == 1);
      BOOST_TEST(pMesh->nVertices() == 2);
      BOOST_TEST(pMesh->vertex(0).getGlobalIndex() == 0);
      BOOST_TEST(pMesh->vertex(1).getGlobalIndex() == 1);
      BOOST_TEST(pMesh->vertex(0).isOwner() == true);
      BOOST_TEST(pMesh->vertex(1).isOwner() == false);
    } else if (context.isRank(1)) { //Secondary rank 2
      BOOST_TEST(pMesh->nVertices() == 1);
      BOOST_TEST(pMesh->vertex(0).getGlobalIndex() == 1);
      BOOST_TEST(pMesh->vertex(0).isOwner() == true);
    } else if (context.isRank(2)) { //Secondary rank 3
      BOOST_TEST(pMesh->nVertices() == 0);
    }
  }
}

BOOST_AUTO_TEST_CASE(ProvideAndReceiveCouplingMode)
{
  PRECICE_TEST("Fluid"_on(1_rank), "Solid"_on(1_rank), Require::Events);
  auto m2n = context.connectPrimaryRanks("Solid", "Fluid");

  int dimensions = 2;

  if (context.isNamed("Solid")) {
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, testing::nextMeshID()));
    createSolidzMesh2D(pSolidzMesh);
    ProvidedPartition part(pSolidzMesh);
    part.addM2N(m2n);
    part.communicate();
    part.compute();

    BOOST_TEST(pSolidzMesh->getGlobalNumberOfVertices() == 6);
    BOOST_TEST(pSolidzMesh->nVertices() == 6);
    BOOST_TEST(pSolidzMesh->edges().size() == 5);
    BOOST_TEST(pSolidzMesh->getVertexOffsets().size() == 1);
    BOOST_TEST(pSolidzMesh->getVertexOffsets().at(0) == 6);
    BOOST_TEST(pSolidzMesh->vertex(0).getGlobalIndex() == 0);
    BOOST_TEST(pSolidzMesh->vertex(1).getGlobalIndex() == 1);
    BOOST_TEST(pSolidzMesh->vertex(2).getGlobalIndex() == 2);
    BOOST_TEST(pSolidzMesh->vertex(3).getGlobalIndex() == 3);
    BOOST_TEST(pSolidzMesh->vertex(4).getGlobalIndex() == 4);
    BOOST_TEST(pSolidzMesh->vertex(5).getGlobalIndex() == 5);
    BOOST_TEST(pSolidzMesh->vertex(0).isOwner() == true);
    BOOST_TEST(pSolidzMesh->vertex(1).isOwner() == true);
    BOOST_TEST(pSolidzMesh->vertex(2).isOwner() == true);
    BOOST_TEST(pSolidzMesh->vertex(3).isOwner() == true);
    BOOST_TEST(pSolidzMesh->vertex(4).isOwner() == true);
    BOOST_TEST(pSolidzMesh->vertex(5).isOwner() == true);
  } else {
    BOOST_TEST(context.isNamed("Fluid"));
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, testing::nextMeshID()));

    mesh::PtrMesh       pOtherMesh(new mesh::Mesh("OtherMesh", dimensions, testing::nextMeshID()));
    mapping::PtrMapping boundingFromMapping = mapping::PtrMapping(
        new mapping::NearestNeighborMapping(mapping::Mapping::CONSISTENT, dimensions));
    boundingFromMapping->setMeshes(pSolidzMesh, pOtherMesh);

    double            safetyFactor = 0.1;
    ReceivedPartition part(pSolidzMesh, ReceivedPartition::ON_PRIMARY_RANK, safetyFactor);
    part.addFromMapping(boundingFromMapping);
    part.addM2N(m2n);
    part.communicate();
    part.compute();

    BOOST_TEST(pSolidzMesh->getGlobalNumberOfVertices() == 6);
    BOOST_TEST(pSolidzMesh->nVertices() == 6);
    BOOST_TEST(pSolidzMesh->edges().size() == 5);
    BOOST_TEST(pSolidzMesh->getVertexOffsets().size() == 1);
    BOOST_TEST(pSolidzMesh->getVertexOffsets().at(0) == 6);
    BOOST_TEST(pSolidzMesh->vertex(0).getGlobalIndex() == 0);
    BOOST_TEST(pSolidzMesh->vertex(1).getGlobalIndex() == 1);
    BOOST_TEST(pSolidzMesh->vertex(2).getGlobalIndex() == 2);
    BOOST_TEST(pSolidzMesh->vertex(3).getGlobalIndex() == 3);
    BOOST_TEST(pSolidzMesh->vertex(4).getGlobalIndex() == 4);
    BOOST_TEST(pSolidzMesh->vertex(5).getGlobalIndex() == 5);
    BOOST_TEST(pSolidzMesh->vertex(0).isOwner() == true);
    BOOST_TEST(pSolidzMesh->vertex(1).isOwner() == true);
    BOOST_TEST(pSolidzMesh->vertex(2).isOwner() == true);
    BOOST_TEST(pSolidzMesh->vertex(3).isOwner() == true);
    BOOST_TEST(pSolidzMesh->vertex(4).isOwner() == true);
    BOOST_TEST(pSolidzMesh->vertex(5).isOwner() == true);
  }
}

BOOST_AUTO_TEST_CASE(TestCompareBoundingBoxes2D)
{
  PRECICE_TEST("SOLIDZ"_on(1_rank), "NASTIN"_on(3_ranks).setupIntraComm(), Require::Events);

  testing::ConnectionOptions options;
  options.useOnlyPrimaryCom = false;
  options.useTwoLevelInit   = true;
  auto m2n                  = context.connectPrimaryRanks("SOLIDZ", "NASTIN", options);

  int dimensions = 2;

  // construct send global boundingbox
  mesh::Mesh::BoundingBoxMap sendGlobalBB;
  for (int remoteRank = 0; remoteRank < 3; remoteRank++) {
    std::vector<double> bounds;
    for (int i = 0; i < dimensions; i++) {
      bounds.push_back(3 - remoteRank - 1);
      bounds.push_back(3 - remoteRank);
    }
    sendGlobalBB.emplace(remoteRank, mesh::BoundingBox(bounds));
  }

  if (context.isNamed("SOLIDZ")) {
    int                             connectionMapSize = 0;
    std::map<int, std::vector<int>> receivedConnectionMap;
    mesh::PtrMesh                   pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, testing::nextMeshID()));
    m2n->getPrimaryRankCommunication()->send(3, 0);
    com::sendBoundingBoxMap(*m2n->getPrimaryRankCommunication(), 0, sendGlobalBB);
    std::vector<int> connectedRanksList = m2n->getPrimaryRankCommunication()->receiveRange(0, com::AsVectorTag<int>{});
    connectionMapSize                   = connectedRanksList.size();
    BOOST_TEST_REQUIRE(connectionMapSize == 2);

    std::vector<int> connectedRanks;
    connectedRanks.push_back(-1);
    for (auto &rank : connectedRanksList) {
      receivedConnectionMap[rank] = connectedRanks;
    }

    com::receiveConnectionMap(*m2n->getPrimaryRankCommunication(), 0, receivedConnectionMap);

    // test whether we receive correct connection map
    BOOST_TEST(receivedConnectionMap.at(0).at(0) == 2);
    BOOST_TEST(receivedConnectionMap.at(2).at(0) == 0);

  } else {
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, testing::nextMeshID()));
    mesh::PtrMesh pNastinMesh(new mesh::Mesh("SolidzMesh", dimensions, testing::nextMeshID()));

    mapping::PtrMapping boundingFromMapping = mapping::PtrMapping(
        new mapping::NearestNeighborMapping(mapping::Mapping::CONSISTENT, dimensions));
    mapping::PtrMapping boundingToMapping = mapping::PtrMapping(
        new mapping::NearestNeighborMapping(mapping::Mapping::CONSERVATIVE, dimensions));

    boundingFromMapping->setMeshes(pSolidzMesh, pNastinMesh);
    boundingToMapping->setMeshes(pNastinMesh, pSolidzMesh);

    createNastinMesh2D2(pNastinMesh, context.rank);

    double safetyFactor = 0.0;

    ReceivedPartition part(pSolidzMesh, ReceivedPartition::NO_FILTER, safetyFactor);
    part.addM2N(m2n);
    part.addFromMapping(boundingFromMapping);
    part.addToMapping(boundingToMapping);
    part.compareBoundingBoxes();
  }
}

BOOST_AUTO_TEST_CASE(TestCompareBoundingBoxes3D)
{
  PRECICE_TEST("SOLIDZ"_on(1_rank), "NASTIN"_on(3_ranks).setupIntraComm(), Require::Events);

  testing::ConnectionOptions options;
  options.useOnlyPrimaryCom = false;
  options.useTwoLevelInit   = true;
  auto m2n                  = context.connectPrimaryRanks("SOLIDZ", "NASTIN", options);

  int dimensions = 3;

  // construct send global boundingbox
  mesh::Mesh::BoundingBoxMap sendGlobalBB;
  for (int remoteRank = 0; remoteRank < 3; remoteRank++) {
    std::vector<double> bounds;
    for (int i = 0; i < dimensions; i++) {
      bounds.push_back(3 - remoteRank - 1);
      bounds.push_back(3 - remoteRank);
    }
    sendGlobalBB.emplace(remoteRank, mesh::BoundingBox(bounds));
  }

  if (context.isNamed("SOLIDZ")) {
    int                             connectionMapSize = 0;
    std::map<int, std::vector<int>> receivedConnectionMap;
    mesh::PtrMesh                   pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, testing::nextMeshID()));
    m2n->getPrimaryRankCommunication()->send(3, 0);
    com::sendBoundingBoxMap(*m2n->getPrimaryRankCommunication(), 0, sendGlobalBB);
    std::vector<int> connectedRanksList = m2n->getPrimaryRankCommunication()->receiveRange(0, com::AsVectorTag<int>{});
    connectionMapSize                   = connectedRanksList.size();
    BOOST_TEST(connectionMapSize == 2);

    std::vector<int> connectedRanks;
    connectedRanks.push_back(-1);
    for (auto &rank : connectedRanksList) {
      receivedConnectionMap[rank] = connectedRanks;
    }

    com::receiveConnectionMap(*m2n->getPrimaryRankCommunication(), 0, receivedConnectionMap);

    // test whether we receive correct connection map
    BOOST_TEST(receivedConnectionMap.at(0).at(0) == 2);
    BOOST_TEST(receivedConnectionMap.at(2).at(0) == 0);

  } else {
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, testing::nextMeshID()));
    mesh::PtrMesh pNastinMesh(new mesh::Mesh("SolidzMesh", dimensions, testing::nextMeshID()));

    mapping::PtrMapping boundingFromMapping = mapping::PtrMapping(
        new mapping::NearestNeighborMapping(mapping::Mapping::CONSISTENT, dimensions));
    mapping::PtrMapping boundingToMapping = mapping::PtrMapping(
        new mapping::NearestNeighborMapping(mapping::Mapping::CONSERVATIVE, dimensions));

    boundingFromMapping->setMeshes(pSolidzMesh, pNastinMesh);
    boundingToMapping->setMeshes(pNastinMesh, pSolidzMesh);

    createNastinMesh3D2(pNastinMesh, context.rank);

    double safetyFactor = 0.0;

    ReceivedPartition part(pSolidzMesh, ReceivedPartition::NO_FILTER, safetyFactor);
    part.addM2N(m2n);
    part.addFromMapping(boundingFromMapping);
    part.addToMapping(boundingToMapping);
    part.compareBoundingBoxes();
  }
}

void testParallelSetOwnerInformation(mesh::PtrMesh mesh, int dimensions)
{
  double safetyFactor = 0;

  testing::ConnectionOptions options;
  options.useOnlyPrimaryCom = false;
  options.useTwoLevelInit   = true;
  options.type              = testing::ConnectionType::PointToPoint;

  auto                                      participantCom = com::PtrCommunication(new com::SocketCommunication());
  m2n::DistributedComFactory::SharedPointer distrFactory;

  auto m2n = m2n::PtrM2N(new m2n::M2N(participantCom, distrFactory, options.useOnlyPrimaryCom, options.useTwoLevelInit));

  mapping::PtrMapping boundingFromMapping = mapping::PtrMapping(new mapping::NearestNeighborMapping(mapping::Mapping::CONSISTENT, dimensions));
  mapping::PtrMapping boundingToMapping   = mapping::PtrMapping(new mapping::NearestNeighborMapping(mapping::Mapping::CONSERVATIVE, dimensions));
  boundingFromMapping->setMeshes(mesh, mesh);
  boundingToMapping->setMeshes(mesh, mesh);

  ReceivedPartition part(mesh, ReceivedPartition::ON_SECONDARY_RANKS, safetyFactor);
  part.addM2N(m2n);

  part.addFromMapping(boundingFromMapping);
  part.addToMapping(boundingToMapping);

  mesh->computeBoundingBox();

  using Access = ReceivedPartitionFixture;
  Access::prepareBoundingBox(part);
  Access::tagMeshFirstRound(part);
  Access::createOwnerInformation(part);
}

BOOST_AUTO_TEST_CASE(parallelSetOwnerInformationVertexCount)
{
  /*
    This test examines an edge case for parallel setOwnerinformation function in receivedpartition.cpp
    for 2LI. The provided mesh includes a vertex at point (0, 0). Initially, all receiving ranks receive
    this vertex, but only one of them can own it. Since the rank 2, has the lowest number of vertices,
    this vertex must belong only to it finally.
  */
  PRECICE_TEST(""_on(4_ranks).setupIntraComm(), Require::Events);
  //mesh creation
  int           dimensions = 2;
  mesh::PtrMesh mesh(new mesh::Mesh("mesh", dimensions, testing::nextMeshID()));

  if (context.isRank(0)) {
    Eigen::VectorXd position(dimensions);
    position << 0.0, 0.0;
    mesh->createVertex(position);
    position << 1.0, 0.0;
    mesh->createVertex(position);
    position << 2.0, 0.0;
    mesh->createVertex(position);
    position << 0.0, 1.0;
    mesh->createVertex(position);
    position << 1.0, 1.0;
    mesh->createVertex(position);
    position << 2.0, 1.0;
    mesh->createVertex(position);

  } else if (context.isRank(1)) {
    Eigen::VectorXd position(dimensions);
    position << 0.0, 0.0;
    mesh->createVertex(position);
    position << -1.0, 0.0;
    mesh->createVertex(position);
    position << -2.0, 0.0;
    mesh->createVertex(position);
    position << -0.5, 1.0;
    mesh->createVertex(position);
    position << -1.0, 1.0;
    mesh->createVertex(position);
    position << -2.0, 1.0;
    mesh->createVertex(position);
  } else if (context.isRank(2)) {
    Eigen::VectorXd position(dimensions);
    position << 0.0, 0.0;
    mesh->createVertex(position);
    position << -1.0, -0.5;
    mesh->createVertex(position);
  } else {
    Eigen::VectorXd position(dimensions);
    position << 0.0, 0.0;
    mesh->createVertex(position);
    position << 1.0, -0.5;
    mesh->createVertex(position);
    position << 2.0, -0.5;
    mesh->createVertex(position);
    position << 0.5, -1.0;
    mesh->createVertex(position);
    position << 1.0, -1.0;
    mesh->createVertex(position);
    position << 2.0, -1.0;
    mesh->createVertex(position);
  }

  mesh->computeBoundingBox();
  mesh->setGlobalNumberOfVertices(mesh->nVertices());

  for (auto &vertex : mesh->vertices()) {
    vertex.setGlobalIndex(vertex.getID() + 5 * utils::IntraComm::getRank());

    if (vertex.coord(0) == 0 && vertex.coord(1) == 0) {
      vertex.setGlobalIndex(0);
    }
  }

  testParallelSetOwnerInformation(mesh, dimensions);

  // to check if all ranks have received the vertex at (0, 0)
  bool includeVertex = false;

  for (auto &vertex : mesh->vertices()) {
    if (vertex.getGlobalIndex() == 0) {
      includeVertex = true;
      if (context.isRank(2)) {
        BOOST_TEST(vertex.isOwner() == 1);
      } else {
        BOOST_TEST(vertex.isOwner() == 0);
      }
    }
    BOOST_TEST(includeVertex == true);
  }
}

BOOST_AUTO_TEST_CASE(parallelSetOwnerInformationLowerRank)
{
  /*
    This test examines an edge case for parallel setOwnerinformation function in receivedpartition.cpp
    for 2LI. The provided mesh includes a vertices at point (0, 0, 0) and (0, 0, 1). Initially, all
    receiving ranks receive this vertex, but only one of them can own it. Since the rank 0, has the lowest
    rank number, this vertex must belong only to this rank.
   */
  PRECICE_TEST(""_on(4_ranks).setupIntraComm(), Require::Events);
  //mesh creation
  int           dimensions = 3;
  mesh::PtrMesh mesh(new mesh::Mesh("mesh", dimensions, testing::nextMeshID()));

  if (context.isRank(0)) {
    Eigen::VectorXd position(dimensions);
    position << 0.0, 0.0, 0.0;
    mesh->createVertex(position);
    position << 1.0, 0.0, 0.0;
    mesh->createVertex(position);
    position << 2.0, 0.0, 0.0;
    mesh->createVertex(position);
    position << 0.0, 1.0, 0.0;
    mesh->createVertex(position);
    position << 1.0, 1.0, 0.0;
    mesh->createVertex(position);
    position << 2.0, 1.0, 0.0;
    mesh->createVertex(position);
    position << 0.0, 0.0, 1.0;
    mesh->createVertex(position);
    position << 1.0, 0.0, 1.0;
    mesh->createVertex(position);
    position << 2.0, 0.0, 1.0;
    mesh->createVertex(position);
    position << 0.0, 1.0, 1.0;
    mesh->createVertex(position);
    position << 1.0, 1.0, 1.0;
    mesh->createVertex(position);
    position << 2.0, 1.0, 1.0;
    mesh->createVertex(position);

  } else if (context.isRank(1)) {
    Eigen::VectorXd position(dimensions);
    position << 0.0, 0.0, 0.0;
    mesh->createVertex(position);
    position << -1.0, 0.0, 0.0;
    mesh->createVertex(position);
    position << -2.0, 0.0, 0.0;
    mesh->createVertex(position);
    position << -0.1, 1.0, 0.0;
    mesh->createVertex(position);
    position << -1.0, 1.0, 0.0;
    mesh->createVertex(position);
    position << -2.0, 1.0, 0.0;
    mesh->createVertex(position);
    position << 0.0, 0.0, 1.0;
    mesh->createVertex(position);
    position << -1.0, 0.0, 1.0;
    mesh->createVertex(position);
    position << -2.0, 0.0, 1.0;
    mesh->createVertex(position);
    position << -0.1, 1.0, 1.0;
    mesh->createVertex(position);
    position << -1.0, 1.0, 1.0;
    mesh->createVertex(position);
    position << -2.0, 1.0, 1.0;
    mesh->createVertex(position);
  } else if (context.isRank(2)) {
    Eigen::VectorXd position(dimensions);
    position << 0.0, 0.0, 0.0;
    mesh->createVertex(position);
    position << -1.0, -0.1, 0.0;
    mesh->createVertex(position);
    position << -2.0, -0.1, 0.0;
    mesh->createVertex(position);
    position << 0.0, -1.0, 0.0;
    mesh->createVertex(position);
    position << -1.0, -1.0, 0.0;
    mesh->createVertex(position);
    position << -2.0, -1.0, 0.0;
    mesh->createVertex(position);
    position << 0.0, 0.0, 1.0;
    mesh->createVertex(position);
    position << -1.0, -0.1, 1.0;
    mesh->createVertex(position);
    position << -2.0, -0.1, 1.0;
    mesh->createVertex(position);
    position << 0.0, -1.0, 1.0;
    mesh->createVertex(position);
    position << -1.0, -1.0, 1.0;
    mesh->createVertex(position);
    position << -2.0, -1.0, 1.0;
    mesh->createVertex(position);
  } else {
    Eigen::VectorXd position(dimensions);
    position << 0.0, 0.0, 0.0;
    mesh->createVertex(position);
    position << 1.0, -0.1, 0.0;
    mesh->createVertex(position);
    position << 2.0, -0.1, 0.0;
    mesh->createVertex(position);
    position << 0.0, -1.0, 0.0;
    mesh->createVertex(position);
    position << 1.0, -1.0, 0.0;
    mesh->createVertex(position);
    position << 2.0, -1.0, 0.0;
    mesh->createVertex(position);
    position << 0.0, 0.0, 1.0;
    mesh->createVertex(position);
    position << 1.0, -0.1, 1.0;
    mesh->createVertex(position);
    position << 2.0, -0.1, 1.0;
    mesh->createVertex(position);
    position << 0.0, -1.0, 1.0;
    mesh->createVertex(position);
    position << 1.0, -1.0, 1.0;
    mesh->createVertex(position);
    position << 2.0, -1.0, 1.0;
    mesh->createVertex(position);
  }

  mesh->computeBoundingBox();
  mesh->setGlobalNumberOfVertices(mesh->nVertices());

  for (auto &vertex : mesh->vertices()) {
    vertex.setGlobalIndex(vertex.getID() + 10 * utils::IntraComm::getRank());

    if (vertex.coord(0) == 0 && vertex.coord(1) == 0) {
      if (vertex.coord(2) == 0) {
        vertex.setGlobalIndex(0);
      } else if (vertex.coord(2) == 1) {
        vertex.setGlobalIndex(6);
      }
    }
  }

  testParallelSetOwnerInformation(mesh, dimensions);

  // to check if all ranks have received the vertex at (0, 0, 0)
  bool includeVertex = false;

  for (auto &vertex : mesh->vertices()) {
    if (vertex.getGlobalIndex() == 0) {
      includeVertex = true;
      if (context.isRank(0)) {
        BOOST_TEST(vertex.isOwner() == 1);
      } else {
        BOOST_TEST(vertex.isOwner() == 0);
      }
    }
    BOOST_TEST(includeVertex == true);
  }
}

BOOST_AUTO_TEST_CASE(parallelSetOwnerInformationEmptyPartition)
{
  /*
    This test examines an edge case for parallel setOwnerinformation function in receivedpartition.cpp
    for 2LI. The provided mesh includes vertices at points (0, 0, 0) and (0, 0, 1). Rank 2 has an
    empty mesh partition. Initially, all ranks (except rank 2) receive this vertex, but only one of them
    can own it. Since the rank 0, has the lowest rank number, this vertex must belong only to this rank.
   */
  PRECICE_TEST(""_on(4_ranks).setupIntraComm(), Require::Events);
  //mesh creation
  int           dimensions = 3;
  mesh::PtrMesh mesh(new mesh::Mesh("mesh", dimensions, testing::nextMeshID()));

  if (context.isRank(0)) {
    Eigen::VectorXd position(dimensions);
    position << 0.0, 0.0, 0.0;
    mesh->createVertex(position);
    position << 1.0, 0.0, 0.0;
    mesh->createVertex(position);
    position << 2.0, 0.0, 0.0;
    mesh->createVertex(position);
    position << 0.0, 1.0, 0.0;
    mesh->createVertex(position);
    position << 1.0, 1.0, 0.0;
    mesh->createVertex(position);
    position << 2.0, 1.0, 0.0;
    mesh->createVertex(position);
    position << 0.0, 0.0, 1.0;
    mesh->createVertex(position);
    position << 1.0, 0.0, 1.0;
    mesh->createVertex(position);
    position << 2.0, 0.0, 1.0;
    mesh->createVertex(position);
    position << 0.0, 1.0, 1.0;
    mesh->createVertex(position);
    position << 1.0, 1.0, 1.0;
    mesh->createVertex(position);
    position << 2.0, 1.0, 1.0;
    mesh->createVertex(position);

  } else if (context.isRank(1)) {
    Eigen::VectorXd position(dimensions);
    position << 0.0, 0.0, 0.0;
    mesh->createVertex(position);
    position << -1.0, 0.0, 0.0;
    mesh->createVertex(position);
    position << -2.0, 0.0, 0.0;
    mesh->createVertex(position);
    position << -0.1, 1.0, 0.0;
    mesh->createVertex(position);
    position << -1.0, 1.0, 0.0;
    mesh->createVertex(position);
    position << -2.0, 1.0, 0.0;
    mesh->createVertex(position);
    position << 0.0, 0.0, 1.0;
    mesh->createVertex(position);
    position << -1.0, 0.0, 1.0;
    mesh->createVertex(position);
    position << -2.0, 0.0, 1.0;
    mesh->createVertex(position);
    position << -0.1, 1.0, 1.0;
    mesh->createVertex(position);
    position << -1.0, 1.0, 1.0;
    mesh->createVertex(position);
    position << -2.0, 1.0, 1.0;
    mesh->createVertex(position);
  } else if (context.isRank(2)) {
  } else {
    Eigen::VectorXd position(dimensions);
    position << 0.0, 0.0, 0.0;
    mesh->createVertex(position);
    position << 1.0, -0.1, 0.0;
    mesh->createVertex(position);
    position << 2.0, -0.1, 0.0;
    mesh->createVertex(position);
    position << 0.0, -1.0, 0.0;
    mesh->createVertex(position);
    position << 1.0, -1.0, 0.0;
    mesh->createVertex(position);
    position << 2.0, -1.0, 0.0;
    mesh->createVertex(position);
    position << 0.0, 0.0, 1.0;
    mesh->createVertex(position);
    position << 1.0, -0.1, 1.0;
    mesh->createVertex(position);
    position << 2.0, -0.1, 1.0;
    mesh->createVertex(position);
    position << 0.0, -1.0, 1.0;
    mesh->createVertex(position);
    position << 1.0, -1.0, 1.0;
    mesh->createVertex(position);
    position << 2.0, -1.0, 1.0;
    mesh->createVertex(position);
  }

  mesh->computeBoundingBox();
  mesh->setGlobalNumberOfVertices(mesh->nVertices());

  for (auto &vertex : mesh->vertices()) {
    vertex.setGlobalIndex(vertex.getID() + 10 * utils::IntraComm::getRank());

    if (vertex.coord(0) == 0 && vertex.coord(1) == 0) {
      if (vertex.coord(2) == 0) {
        vertex.setGlobalIndex(0);
      } else if (vertex.coord(2) == 1) {
        vertex.setGlobalIndex(6);
      }
    }
  }

  testParallelSetOwnerInformation(mesh, dimensions);

  // to check if all ranks have received the vertex at (0, 0, 0)
  bool includeVertex = false;

  for (auto &vertex : mesh->vertices()) {
    if (vertex.getGlobalIndex() == 0) {
      includeVertex = true;
      if (context.isRank(0)) {
        BOOST_TEST(vertex.isOwner() == 1);
      } else {
        BOOST_TEST(vertex.isOwner() == 0);
      }
    }
    BOOST_TEST(includeVertex == true);
  }
}

// Test with two "from" and two "to" mappings
BOOST_AUTO_TEST_CASE(RePartitionMultipleMappings)
{
  PRECICE_TEST("Solid"_on(1_rank), "Fluid"_on(3_ranks).setupIntraComm(), Require::Events);
  auto m2n = context.connectPrimaryRanks("Solid", "Fluid");

  int             dimensions = 2;
  Eigen::VectorXd offset     = Eigen::VectorXd::Zero(dimensions);

  if (context.isNamed("Solid")) { //SOLIDZ
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, testing::nextMeshID()));
    createSolidzMesh2D(pSolidzMesh);
    BOOST_TEST(pSolidzMesh->nVertices() == 6);
    ProvidedPartition part(pSolidzMesh);
    part.addM2N(m2n);
    part.communicate();
  } else {
    BOOST_TEST(context.isNamed("Fluid"));
    mesh::PtrMesh pNastinMesh1(new mesh::Mesh("NastinMesh1", dimensions, testing::nextMeshID()));
    mesh::PtrMesh pNastinMesh2(new mesh::Mesh("NastinMesh1", dimensions, testing::nextMeshID()));
    mesh::PtrMesh pNastinMesh3(new mesh::Mesh("NastinMesh1", dimensions, testing::nextMeshID()));
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, testing::nextMeshID()));

    mapping::PtrMapping boundingFromMapping1 = mapping::PtrMapping(
        new mapping::NearestNeighborMapping(mapping::Mapping::CONSISTENT, dimensions));
    mapping::PtrMapping boundingToMapping1 = mapping::PtrMapping(
        new mapping::NearestNeighborMapping(mapping::Mapping::CONSERVATIVE, dimensions));
    mapping::PtrMapping boundingFromMapping2 = mapping::PtrMapping(
        new mapping::NearestNeighborMapping(mapping::Mapping::CONSISTENT, dimensions));
    mapping::PtrMapping boundingToMapping2 = mapping::PtrMapping(
        new mapping::NearestNeighborMapping(mapping::Mapping::CONSERVATIVE, dimensions));
    boundingFromMapping1->setMeshes(pSolidzMesh, pNastinMesh1);
    boundingToMapping1->setMeshes(pNastinMesh1, pSolidzMesh);
    boundingFromMapping2->setMeshes(pSolidzMesh, pNastinMesh2);
    boundingToMapping2->setMeshes(pNastinMesh3, pSolidzMesh);

    if (context.rank == 0) {
      Eigen::VectorXd position(dimensions);
      position << 0.0, 0.0;
      pNastinMesh1->createVertex(position);
      position << 0.0, 2.0;
      pNastinMesh2->createVertex(position);
    } else if (context.rank == 1) {
      // not at interface
    } else if (context.rank == 2) {
      Eigen::VectorXd position(dimensions);
      position << 0.0, 4.0;
      pNastinMesh2->createVertex(position);
      position << 0.0, 6.0;
      pNastinMesh3->createVertex(position);
    }
    pNastinMesh1->computeBoundingBox();
    pNastinMesh2->computeBoundingBox();
    pNastinMesh3->computeBoundingBox();

    double safetyFactor = 0.1;

    ReceivedPartition part(pSolidzMesh, ReceivedPartition::ON_SECONDARY_RANKS, safetyFactor);
    part.addM2N(m2n);
    part.addFromMapping(boundingFromMapping1);
    part.addToMapping(boundingToMapping1);
    part.addFromMapping(boundingFromMapping2);
    part.addToMapping(boundingToMapping2);
    part.communicate();
    part.compute();

    BOOST_TEST_CONTEXT(*pSolidzMesh)
    {
      // check if the sending and filtering worked right
      if (context.isPrimary()) { //Primary
        BOOST_TEST(pSolidzMesh->nVertices() == 2);
        BOOST_TEST(pSolidzMesh->edges().size() == 1);
      } else if (context.isRank(1)) { //SecondaryRank1
        BOOST_TEST(pSolidzMesh->nVertices() == 0);
        BOOST_TEST(pSolidzMesh->edges().size() == 0);
      } else if (context.isRank(2)) { //Secondary rank 2
        BOOST_TEST(pSolidzMesh->nVertices() == 2);
        BOOST_TEST(pSolidzMesh->edges().size() == 1);
      }
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

#endif // PRECICE_NO_MPI
