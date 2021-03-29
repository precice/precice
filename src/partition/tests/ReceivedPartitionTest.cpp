#ifndef PRECICE_NO_MPI

#include <Eigen/Core>
#include <algorithm>
#include <map>
#include <memory>
#include <string>
#include <vector>
#include "com/CommunicateBoundingBox.hpp"
#include "com/Communication.hpp"
#include "com/SharedPointer.hpp"
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
#include "precice/impl/versions.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"
#include "utils/assertion.hpp"

namespace precice {
namespace mesh {
class Edge;
} // namespace mesh
} // namespace precice

using namespace precice;
using namespace partition;
using precice::testing::TestContext;

BOOST_AUTO_TEST_SUITE(PartitionTests)
BOOST_AUTO_TEST_SUITE(ReceivedPartitionTests)

void tearDownParallelEnvironment()
{
  mesh::Data::resetDataCount();
}

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
  pSolidzMesh->computeState();
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
  pSolidzMesh->computeState();
  pSolidzMesh->computeBoundingBox();
}

void createNastinMesh2D(mesh::PtrMesh pNastinMesh, int rank)
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
  pNastinMesh->computeState();
  pNastinMesh->computeBoundingBox();
}

void createNastinMesh2D2(mesh::PtrMesh pNastinMesh, int rank)
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
  pNastinMesh->computeState();
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
  pSolidzMesh->computeState();
  pSolidzMesh->computeBoundingBox();
}

void createNastinMesh3D(mesh::PtrMesh pNastinMesh, int rank)
{
  int dimensions = 3;
  BOOST_TEST(pNastinMesh);
  BOOST_TEST(pNastinMesh->getDimensions() == dimensions);

  if (rank == 0) { //Master
    Eigen::VectorXd position(dimensions);
    position << -1.0, -1.0, 0.0;
    pNastinMesh->createVertex(position);
    position << -0.75, -0.75, 0.5;
    pNastinMesh->createVertex(position);
  } else if (rank == 1) { //Slave1
    // slave1 not at interface
  } else if (rank == 2) { //Slave2
    Eigen::VectorXd position(dimensions);
    position << 0.0, 0.0, -1.0;
    pNastinMesh->createVertex(position);
    position << 0.5, 0.5, 0.0;
    pNastinMesh->createVertex(position);
  }
  pNastinMesh->computeState();
  pNastinMesh->computeBoundingBox();
}

void createNastinMesh3D2(mesh::PtrMesh pNastinMesh, int rank)
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
  pNastinMesh->computeState();
  pNastinMesh->computeBoundingBox();
}

BOOST_AUTO_TEST_CASE(RePartitionNNBroadcastFilter2D)
{
  PRECICE_TEST("Solid"_on(1_rank), "Fluid"_on(3_ranks).setupMasterSlaves(), Require::Events);
  auto m2n = context.connectMasters("Solid", "Fluid");

  int             dimensions  = 2;
  bool            flipNormals = false;
  Eigen::VectorXd offset      = Eigen::VectorXd::Zero(dimensions);

  if (context.isNamed("Solid")) { //SOLIDZ
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals, testing::nextMeshID()));
    createSolidzMesh2D(pSolidzMesh);
    BOOST_TEST(pSolidzMesh->vertices().size() == 6);
    ProvidedPartition part(pSolidzMesh);
    part.addM2N(m2n);
    part.communicate();
  } else {
    BOOST_TEST(context.isNamed("Fluid"));
    mesh::PtrMesh pNastinMesh(new mesh::Mesh("NastinMesh", dimensions, flipNormals, testing::nextMeshID()));
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals, testing::nextMeshID()));

    mapping::PtrMapping boundingFromMapping = mapping::PtrMapping(
        new mapping::NearestNeighborMapping(mapping::Mapping::CONSISTENT, dimensions));
    mapping::PtrMapping boundingToMapping = mapping::PtrMapping(
        new mapping::NearestNeighborMapping(mapping::Mapping::CONSERVATIVE, dimensions));
    boundingFromMapping->setMeshes(pSolidzMesh, pNastinMesh);
    boundingToMapping->setMeshes(pNastinMesh, pSolidzMesh);

    createNastinMesh2D(pNastinMesh, context.rank);

    double safetyFactor = 0.1;

    ReceivedPartition part(pSolidzMesh, ReceivedPartition::ON_MASTER, safetyFactor);
    part.addM2N(m2n);
    part.addFromMapping(boundingFromMapping);
    part.addToMapping(boundingToMapping);
    part.communicate();
    part.compute();

    BOOST_TEST_CONTEXT(*pSolidzMesh)
    {
      // check if the sending and filtering worked right
      if (context.isMaster()) { //Master
        BOOST_TEST(pSolidzMesh->vertices().size() == 2);
        BOOST_TEST(pSolidzMesh->edges().size() == 1);
      } else if (context.isRank(1)) { //Slave1
        BOOST_TEST(pSolidzMesh->vertices().size() == 0);
        BOOST_TEST(pSolidzMesh->edges().size() == 0);
      } else if (context.isRank(2)) { //Slave2
        BOOST_TEST(pSolidzMesh->vertices().size() == 2);
        BOOST_TEST(pSolidzMesh->edges().size() == 1);
      }
    }
  }

  tearDownParallelEnvironment();
}

BOOST_AUTO_TEST_CASE(RePartitionNNDoubleNode2D)
{
  PRECICE_TEST("Solid"_on(1_rank), "Fluid"_on(3_ranks).setupMasterSlaves(), Require::Events);
  auto m2n = context.connectMasters("Solid", "Fluid");

  int             dimensions  = 2;
  bool            flipNormals = false;
  Eigen::VectorXd offset      = Eigen::VectorXd::Zero(dimensions);

  if (context.isNamed("Solid")) { //SOLIDZ
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals, testing::nextMeshID()));
    createSolidzMesh2DSmall(pSolidzMesh);
    ProvidedPartition part(pSolidzMesh);
    part.addM2N(m2n);
    part.communicate();
  } else {
    BOOST_TEST(context.isNamed("Fluid"));
    mesh::PtrMesh pNastinMesh(new mesh::Mesh("NastinMesh", dimensions, flipNormals, testing::nextMeshID()));
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals, testing::nextMeshID()));

    mapping::PtrMapping boundingFromMapping = mapping::PtrMapping(
        new mapping::NearestNeighborMapping(mapping::Mapping::CONSISTENT, dimensions));
    mapping::PtrMapping boundingToMapping = mapping::PtrMapping(
        new mapping::NearestNeighborMapping(mapping::Mapping::CONSERVATIVE, dimensions));
    boundingFromMapping->setMeshes(pSolidzMesh, pNastinMesh);
    boundingToMapping->setMeshes(pNastinMesh, pSolidzMesh);

    createNastinMesh2D(pNastinMesh, context.rank);

    double safetyFactor = 0.5;

    ReceivedPartition part(pSolidzMesh, ReceivedPartition::ON_SLAVES, safetyFactor);
    part.addM2N(m2n);
    part.addFromMapping(boundingFromMapping);
    part.addToMapping(boundingToMapping);
    part.communicate();
    part.compute();

    // check if the sending and filtering worked right
    if (context.isMaster()) { //Master
      BOOST_TEST(pSolidzMesh->vertices().size() == 2);
      BOOST_TEST(pSolidzMesh->edges().size() == 1);
    } else if (context.isRank(1)) { //Slave1
      BOOST_TEST(pSolidzMesh->vertices().size() == 0);
      BOOST_TEST(pSolidzMesh->edges().size() == 0);
    } else if (context.isRank(2)) { //Slave2
      BOOST_TEST(pSolidzMesh->vertices().size() == 2);
      BOOST_TEST(pSolidzMesh->edges().size() == 1);
    }
  }
  tearDownParallelEnvironment();
}

BOOST_AUTO_TEST_CASE(RePartitionNPPreFilterPostFilter2D)
{
  PRECICE_TEST("Solid"_on(1_rank), "Fluid"_on(3_ranks).setupMasterSlaves(), Require::Events);
  auto m2n = context.connectMasters("Solid", "Fluid");

  int  dimensions  = 2;
  bool flipNormals = false;

  if (context.isNamed("Solid")) { //SOLIDZ
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals, testing::nextMeshID()));
    createSolidzMesh2D(pSolidzMesh);
    ProvidedPartition part(pSolidzMesh);
    part.addM2N(m2n);
    part.communicate();
  } else {
    BOOST_TEST(context.isNamed("Fluid"));
    mesh::PtrMesh pNastinMesh(new mesh::Mesh("NastinMesh", dimensions, flipNormals, testing::nextMeshID()));
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals, testing::nextMeshID()));

    mapping::PtrMapping boundingFromMapping = mapping::PtrMapping(
        new mapping::NearestProjectionMapping(mapping::Mapping::CONSISTENT, dimensions));
    mapping::PtrMapping boundingToMapping = mapping::PtrMapping(
        new mapping::NearestNeighborMapping(mapping::Mapping::CONSERVATIVE, dimensions));
    boundingFromMapping->setMeshes(pSolidzMesh, pNastinMesh);
    boundingToMapping->setMeshes(pNastinMesh, pSolidzMesh);

    createNastinMesh2D(pNastinMesh, context.rank);

    double            safetyFactor = 0.1;
    ReceivedPartition part(pSolidzMesh, ReceivedPartition::ON_MASTER, safetyFactor);
    part.addM2N(m2n);
    part.addFromMapping(boundingFromMapping);
    part.addToMapping(boundingToMapping);
    part.communicate();
    part.compute();

    BOOST_TEST_CONTEXT(*pSolidzMesh)
    {
      // check if the sending and filtering worked right
      if (context.isMaster()) { //Master
        BOOST_TEST(pSolidzMesh->vertices().size() == 3);
        BOOST_TEST(pSolidzMesh->edges().size() == 2);
      } else if (context.isRank(1)) { //Slave1
        BOOST_TEST(pSolidzMesh->vertices().size() == 0);
        BOOST_TEST(pSolidzMesh->edges().size() == 0);
      } else if (context.isRank(2)) { //Slave2
        BOOST_TEST(pSolidzMesh->vertices().size() == 3);
        BOOST_TEST(pSolidzMesh->edges().size() == 2);
      }
    }
  }
  tearDownParallelEnvironment();
}

#ifndef PRECICE_NO_PETSC
BOOST_AUTO_TEST_CASE(RePartitionRBFGlobal2D)
{
  PRECICE_TEST("Solid"_on(1_rank), "Fluid"_on(3_ranks).setupMasterSlaves(), Require::Events, Require::PETSc);
  auto m2n = context.connectMasters("Solid", "Fluid");

  int  dimensions  = 2;
  bool flipNormals = false;

  if (context.isNamed("Solid")) { //SOLIDZ
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals, testing::nextMeshID()));
    createSolidzMesh2D(pSolidzMesh);
    ProvidedPartition part(pSolidzMesh);
    part.addM2N(m2n);
    part.communicate();
  } else {
    BOOST_TEST(context.isNamed("Fluid"));
    mesh::PtrMesh pNastinMesh(new mesh::Mesh("NastinMesh", dimensions, flipNormals, testing::nextMeshID()));
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals, testing::nextMeshID()));

    mapping::PtrMapping boundingFromMapping = mapping::PtrMapping(
        new mapping::PetRadialBasisFctMapping<mapping::ThinPlateSplines>(mapping::Mapping::CONSISTENT, dimensions,
                                                                         mapping::ThinPlateSplines(), false, false, false));
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
      if (context.isMaster()) { //Master
        BOOST_TEST(pSolidzMesh->vertices().size() == 6);
        BOOST_TEST(pSolidzMesh->edges().size() == 5);
        BOOST_TEST(pSolidzMesh->vertices().at(0).isOwner() == true);
        BOOST_TEST(pSolidzMesh->vertices().at(1).isOwner() == true);
        BOOST_TEST(pSolidzMesh->vertices().at(2).isOwner() == true);
        BOOST_TEST(pSolidzMesh->vertices().at(3).isOwner() == false);
        BOOST_TEST(pSolidzMesh->vertices().at(4).isOwner() == false);
        BOOST_TEST(pSolidzMesh->vertices().at(5).isOwner() == false);
        BOOST_TEST(pSolidzMesh->vertices().at(0).getGlobalIndex() == 0);
        BOOST_TEST(pSolidzMesh->vertices().at(1).getGlobalIndex() == 1);
        BOOST_TEST(pSolidzMesh->vertices().at(2).getGlobalIndex() == 2);
        BOOST_TEST(pSolidzMesh->vertices().at(3).getGlobalIndex() == 3);
        BOOST_TEST(pSolidzMesh->vertices().at(4).getGlobalIndex() == 4);
        BOOST_TEST(pSolidzMesh->vertices().at(5).getGlobalIndex() == 5);
      } else if (context.isRank(1)) { //Slave2
        BOOST_TEST(pSolidzMesh->vertices().size() == 0);
        BOOST_TEST(pSolidzMesh->edges().size() == 0);
      } else if (context.isRank(2)) { //Slave3
        BOOST_TEST(pSolidzMesh->vertices().size() == 6);
        BOOST_TEST(pSolidzMesh->edges().size() == 5);
        BOOST_TEST(pSolidzMesh->vertices().at(0).isOwner() == false);
        BOOST_TEST(pSolidzMesh->vertices().at(1).isOwner() == false);
        BOOST_TEST(pSolidzMesh->vertices().at(2).isOwner() == false);
        BOOST_TEST(pSolidzMesh->vertices().at(3).isOwner() == true);
        BOOST_TEST(pSolidzMesh->vertices().at(4).isOwner() == true);
        BOOST_TEST(pSolidzMesh->vertices().at(5).isOwner() == true);
        BOOST_TEST(pSolidzMesh->vertices().at(0).getGlobalIndex() == 0);
        BOOST_TEST(pSolidzMesh->vertices().at(1).getGlobalIndex() == 1);
        BOOST_TEST(pSolidzMesh->vertices().at(2).getGlobalIndex() == 2);
        BOOST_TEST(pSolidzMesh->vertices().at(3).getGlobalIndex() == 3);
        BOOST_TEST(pSolidzMesh->vertices().at(4).getGlobalIndex() == 4);
        BOOST_TEST(pSolidzMesh->vertices().at(5).getGlobalIndex() == 5);
      }
    }
  }
  tearDownParallelEnvironment();
}

BOOST_AUTO_TEST_CASE(RePartitionRBFLocal2D1)
{
  PRECICE_TEST("Solid"_on(1_rank), "Fluid"_on(3_ranks).setupMasterSlaves(), Require::Events, Require::PETSc);
  auto m2n = context.connectMasters("Solid", "Fluid");

  int  dimensions  = 2;
  bool flipNormals = false;

  if (context.isNamed("Solid")) { //SOLIDZ
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals, testing::nextMeshID()));
    createSolidzMesh2D(pSolidzMesh);
    ProvidedPartition part(pSolidzMesh);
    part.addM2N(m2n);
    part.communicate();
  } else {
    BOOST_TEST(context.isNamed("Fluid"));
    mesh::PtrMesh pNastinMesh(new mesh::Mesh("NastinMesh", dimensions, flipNormals, testing::nextMeshID()));
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals, testing::nextMeshID()));

    double supportRadius = 0.25;

    mapping::PtrMapping boundingFromMapping = mapping::PtrMapping(
        new mapping::PetRadialBasisFctMapping<mapping::CompactThinPlateSplinesC2>(mapping::Mapping::CONSISTENT, dimensions,
                                                                                  mapping::CompactThinPlateSplinesC2(supportRadius), false, false, false));
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
      if (context.isMaster()) { //Master
        BOOST_TEST(pSolidzMesh->vertices().size() == 3);
        BOOST_TEST(pSolidzMesh->edges().size() == 2);
        BOOST_TEST(pSolidzMesh->vertices().at(0).isOwner() == true);
        BOOST_TEST(pSolidzMesh->vertices().at(1).isOwner() == true);
        BOOST_TEST(pSolidzMesh->vertices().at(2).isOwner() == true);
        BOOST_TEST(pSolidzMesh->vertices().at(0).getGlobalIndex() == 0);
        BOOST_TEST(pSolidzMesh->vertices().at(1).getGlobalIndex() == 1);
        BOOST_TEST(pSolidzMesh->vertices().at(2).getGlobalIndex() == 2);
      } else if (context.isRank(1)) { //Slave2
        BOOST_TEST(pSolidzMesh->vertices().size() == 0);
        BOOST_TEST(pSolidzMesh->edges().size() == 0);
      } else if (context.isRank(2)) { //Slave3
        BOOST_TEST(pSolidzMesh->vertices().size() == 3);
        BOOST_TEST(pSolidzMesh->edges().size() == 2);
        BOOST_TEST(pSolidzMesh->vertices().at(0).isOwner() == true);
        BOOST_TEST(pSolidzMesh->vertices().at(1).isOwner() == true);
        BOOST_TEST(pSolidzMesh->vertices().at(2).isOwner() == true);
        BOOST_TEST(pSolidzMesh->vertices().at(0).getGlobalIndex() == 3);
        BOOST_TEST(pSolidzMesh->vertices().at(1).getGlobalIndex() == 4);
        BOOST_TEST(pSolidzMesh->vertices().at(2).getGlobalIndex() == 5);
      }
    }
  }
  tearDownParallelEnvironment();
}

BOOST_AUTO_TEST_CASE(RePartitionRBFLocal2D2)
{
  PRECICE_TEST("Solid"_on(1_rank), "Fluid"_on(3_ranks).setupMasterSlaves(), Require::Events, Require::PETSc);
  auto m2n = context.connectMasters("Solid", "Fluid");

  int  dimensions  = 2;
  bool flipNormals = false;

  if (context.isNamed("Solid")) { //SOLIDZ
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals, testing::nextMeshID()));
    createSolidzMesh2D(pSolidzMesh);
    ProvidedPartition part(pSolidzMesh);
    part.addM2N(m2n);
    part.communicate();
  } else {
    BOOST_TEST(context.isNamed("Fluid"));
    mesh::PtrMesh pNastinMesh(new mesh::Mesh("NastinMesh", dimensions, flipNormals, testing::nextMeshID()));
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals, testing::nextMeshID()));

    double supportRadius = 2.45;

    mapping::PtrMapping boundingFromMapping = mapping::PtrMapping(
        new mapping::PetRadialBasisFctMapping<mapping::CompactThinPlateSplinesC2>(mapping::Mapping::CONSISTENT, dimensions,
                                                                                  mapping::CompactThinPlateSplinesC2(supportRadius), false, false, false));
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
      if (context.isMaster()) { //Master
        BOOST_TEST(pSolidzMesh->vertices().size() == 4);
        BOOST_TEST(pSolidzMesh->edges().size() == 3);
        BOOST_TEST(pSolidzMesh->vertices().at(0).isOwner() == true);
        BOOST_TEST(pSolidzMesh->vertices().at(1).isOwner() == true);
        BOOST_TEST(pSolidzMesh->vertices().at(2).isOwner() == true);
        BOOST_TEST(pSolidzMesh->vertices().at(3).isOwner() == false);
        BOOST_TEST(pSolidzMesh->vertices().at(0).getGlobalIndex() == 0);
        BOOST_TEST(pSolidzMesh->vertices().at(1).getGlobalIndex() == 1);
        BOOST_TEST(pSolidzMesh->vertices().at(2).getGlobalIndex() == 2);
        BOOST_TEST(pSolidzMesh->vertices().at(3).getGlobalIndex() == 3);
      } else if (context.isRank(1)) { //Slave2
        BOOST_TEST(pSolidzMesh->vertices().size() == 0);
        BOOST_TEST(pSolidzMesh->edges().size() == 0);
      } else if (context.isRank(2)) { //Slave3
        BOOST_TEST(pSolidzMesh->vertices().size() == 5);
        BOOST_TEST(pSolidzMesh->edges().size() == 4);
        BOOST_TEST(pSolidzMesh->vertices().at(0).isOwner() == false);
        BOOST_TEST(pSolidzMesh->vertices().at(1).isOwner() == false);
        BOOST_TEST(pSolidzMesh->vertices().at(2).isOwner() == true);
        BOOST_TEST(pSolidzMesh->vertices().at(3).isOwner() == true);
        BOOST_TEST(pSolidzMesh->vertices().at(4).isOwner() == true);
        BOOST_TEST(pSolidzMesh->vertices().at(0).getGlobalIndex() == 1);
        BOOST_TEST(pSolidzMesh->vertices().at(1).getGlobalIndex() == 2);
        BOOST_TEST(pSolidzMesh->vertices().at(2).getGlobalIndex() == 3);
        BOOST_TEST(pSolidzMesh->vertices().at(3).getGlobalIndex() == 4);
        BOOST_TEST(pSolidzMesh->vertices().at(4).getGlobalIndex() == 5);
      }
    }
  }
  tearDownParallelEnvironment();
}

BOOST_AUTO_TEST_CASE(RePartitionRBFLocal3D)
{
  PRECICE_TEST("Solid"_on(1_rank), "Fluid"_on(3_ranks).setupMasterSlaves(), Require::Events, Require::PETSc);
  auto m2n = context.connectMasters("Solid", "Fluid");

  int  dimensions  = 3;
  bool flipNormals = false;

  if (context.isNamed("Solid")) { //SOLIDZ
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals, testing::nextMeshID()));
    createSolidzMesh3D(pSolidzMesh);
    ProvidedPartition part(pSolidzMesh);
    part.addM2N(m2n);
    part.communicate();
  } else {
    BOOST_TEST(context.isNamed("Fluid"));
    mesh::PtrMesh pNastinMesh(new mesh::Mesh("NastinMesh", dimensions, flipNormals, testing::nextMeshID()));
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals, testing::nextMeshID()));

    double supportRadius1 = 1.2;
    double supportRadius2 = 0.2;

    mapping::PtrMapping boundingFromMapping = mapping::PtrMapping(
        new mapping::PetRadialBasisFctMapping<mapping::CompactThinPlateSplinesC2>(mapping::Mapping::CONSISTENT, dimensions,
                                                                                  mapping::CompactThinPlateSplinesC2(supportRadius1), false, false, false));
    mapping::PtrMapping boundingToMapping = mapping::PtrMapping(
        new mapping::PetRadialBasisFctMapping<mapping::CompactThinPlateSplinesC2>(mapping::Mapping::CONSERVATIVE, dimensions,
                                                                                  mapping::CompactThinPlateSplinesC2(supportRadius2), false, false, false));
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
      if (context.isMaster()) { //Master
        BOOST_TEST(pSolidzMesh->vertices().size() == 5);
        BOOST_TEST(pSolidzMesh->edges().size() == 6);
        BOOST_TEST(pSolidzMesh->triangles().size() == 2);
        BOOST_TEST(pSolidzMesh->vertices().at(0).isOwner() == true);
        BOOST_TEST(pSolidzMesh->vertices().at(1).isOwner() == true);
        BOOST_TEST(pSolidzMesh->vertices().at(2).isOwner() == false);
        BOOST_TEST(pSolidzMesh->vertices().at(3).isOwner() == false);
        BOOST_TEST(pSolidzMesh->vertices().at(4).isOwner() == false);
        BOOST_TEST(pSolidzMesh->vertices().at(0).getGlobalIndex() == 0);
        BOOST_TEST(pSolidzMesh->vertices().at(1).getGlobalIndex() == 1);
        BOOST_TEST(pSolidzMesh->vertices().at(2).getGlobalIndex() == 2);
        BOOST_TEST(pSolidzMesh->vertices().at(3).getGlobalIndex() == 3);
        BOOST_TEST(pSolidzMesh->vertices().at(4).getGlobalIndex() == 4);
      } else if (context.isRank(1)) { //Slave2
        BOOST_TEST(pSolidzMesh->vertices().size() == 0);
        BOOST_TEST(pSolidzMesh->edges().size() == 0);
        BOOST_TEST(pSolidzMesh->triangles().size() == 0);
      } else if (context.isRank(2)) { //Slave3
        BOOST_TEST(pSolidzMesh->vertices().size() == 5);
        BOOST_TEST(pSolidzMesh->edges().size() == 6);
        BOOST_TEST(pSolidzMesh->triangles().size() == 2);
        BOOST_TEST(pSolidzMesh->vertices().at(0).isOwner() == false);
        BOOST_TEST(pSolidzMesh->vertices().at(1).isOwner() == false);
        BOOST_TEST(pSolidzMesh->vertices().at(2).isOwner() == true);
        BOOST_TEST(pSolidzMesh->vertices().at(3).isOwner() == true);
        BOOST_TEST(pSolidzMesh->vertices().at(4).isOwner() == true);
        BOOST_TEST(pSolidzMesh->vertices().at(0).getGlobalIndex() == 0);
        BOOST_TEST(pSolidzMesh->vertices().at(1).getGlobalIndex() == 1);
        BOOST_TEST(pSolidzMesh->vertices().at(2).getGlobalIndex() == 2);
        BOOST_TEST(pSolidzMesh->vertices().at(3).getGlobalIndex() == 3);
        BOOST_TEST(pSolidzMesh->vertices().at(4).getGlobalIndex() == 4);
      }
    }
  }
  tearDownParallelEnvironment();
}

#endif // PRECICE_NO_PETSC

BOOST_AUTO_TEST_CASE(RePartitionNPBroadcastFilter3D)
{
  PRECICE_TEST("Fluid"_on(3_ranks).setupMasterSlaves(), "Solid"_on(1_rank), Require::Events);
  auto m2n = context.connectMasters("Solid", "Fluid");

  int  dimensions  = 3;
  bool flipNormals = false;

  if (context.isNamed("Solid")) { //SOLIDZ
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals, testing::nextMeshID()));
    createSolidzMesh3D(pSolidzMesh);
    ProvidedPartition part(pSolidzMesh);
    part.addM2N(m2n);
    part.communicate();
  } else {
    BOOST_TEST(context.isNamed("Fluid"));
    mesh::PtrMesh pNastinMesh(new mesh::Mesh("NastinMesh", dimensions, flipNormals, testing::nextMeshID()));
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals, testing::nextMeshID()));

    mapping::PtrMapping boundingFromMapping = mapping::PtrMapping(
        new mapping::NearestProjectionMapping(mapping::Mapping::CONSISTENT, dimensions));
    mapping::PtrMapping boundingToMapping = mapping::PtrMapping(
        new mapping::NearestNeighborMapping(mapping::Mapping::CONSERVATIVE, dimensions));
    boundingFromMapping->setMeshes(pSolidzMesh, pNastinMesh);
    boundingToMapping->setMeshes(pNastinMesh, pSolidzMesh);

    createNastinMesh3D(pNastinMesh, context.rank);

    double            safetyFactor = 20.0;
    ReceivedPartition part(pSolidzMesh, ReceivedPartition::ON_MASTER, safetyFactor);
    part.addM2N(m2n);
    part.addFromMapping(boundingFromMapping);
    part.addToMapping(boundingToMapping);
    part.communicate();
    part.compute();

    // check if the sending and filtering worked right
    if (context.isMaster()) { //Master
      BOOST_TEST(pSolidzMesh->vertices().size() == 2);
      BOOST_TEST(pSolidzMesh->edges().size() == 1);
      BOOST_TEST(pSolidzMesh->triangles().size() == 0);
    } else if (context.isRank(1)) { //Slave1
      BOOST_TEST(pSolidzMesh->vertices().size() == 0);
      BOOST_TEST(pSolidzMesh->edges().size() == 0);
      BOOST_TEST(pSolidzMesh->triangles().size() == 0);
    } else if (context.isRank(2)) { //Slave2
      BOOST_TEST(pSolidzMesh->vertices().size() == 3);
      BOOST_TEST(pSolidzMesh->edges().size() == 3);
      BOOST_TEST(pSolidzMesh->triangles().size() == 1);
    }
  }
  tearDownParallelEnvironment();
}

BOOST_AUTO_TEST_CASE(TestRepartitionAndDistribution2D)
{
  PRECICE_TEST("Solid"_on(1_rank), "Fluid"_on(3_ranks).setupMasterSlaves(), Require::Events);
  auto m2n = context.connectMasters("Solid", "Fluid");

  int  dimensions  = 2;
  bool flipNormals = false;

  if (context.isNamed("Solid")) { //SOLIDZ
    mesh::PtrMesh pMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals, testing::nextMeshID()));

    Eigen::VectorXd position(dimensions);
    position << 0.0, 0.0;
    pMesh->createVertex(position);
    position << 1.0, 0.0;
    pMesh->createVertex(position);
    position << 2.0, 0.0;
    pMesh->createVertex(position);

    pMesh->computeState();
    pMesh->computeBoundingBox();

    ProvidedPartition part(pMesh);
    part.addM2N(m2n);
    part.communicate();

  } else {
    BOOST_TEST(context.isNamed("Fluid"));
    mesh::PtrMesh pMesh(new mesh::Mesh("NastinMesh", dimensions, flipNormals, testing::nextMeshID()));
    mesh::PtrMesh pOtherMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals, testing::nextMeshID()));

    mapping::PtrMapping boundingFromMapping = mapping::PtrMapping(
        new mapping::NearestNeighborMapping(mapping::Mapping::CONSISTENT, dimensions));
    boundingFromMapping->setMeshes(pMesh, pOtherMesh);

    if (context.isMaster()) { //Master
      Eigen::VectorXd position(dimensions);
      position << 0.0, 0.0;
      pOtherMesh->createVertex(position);
      position << 0.8, 0.0;
      pOtherMesh->createVertex(position);
    } else if (context.isRank(1)) { //Slave2
      Eigen::VectorXd position(dimensions);
      position << 1.0, 0.0;
      pOtherMesh->createVertex(position);
      position << 1.2, 0.0;
      pOtherMesh->createVertex(position);
    } else if (context.isRank(2)) { //Slave3
      // no vertices
    }

    pOtherMesh->computeState();
    pOtherMesh->computeBoundingBox();

    double            safetyFactor = 20.0;
    ReceivedPartition part(pMesh, ReceivedPartition::ON_MASTER, safetyFactor);
    part.addM2N(m2n);
    part.addFromMapping(boundingFromMapping);
    part.communicate();
    part.compute();

    BOOST_TEST(pMesh->getVertexOffsets().size() == 3);
    BOOST_TEST(pMesh->getVertexOffsets().at(0) == 2);
    BOOST_TEST(pMesh->getVertexOffsets().at(1) == 3);
    BOOST_TEST(pMesh->getVertexOffsets().at(2) == 3);

    if (context.isMaster()) { //Master
      BOOST_TEST(pMesh->getVertexDistribution().at(0).size() == 2);
      BOOST_TEST(pMesh->getVertexDistribution().at(1).size() == 1);
      BOOST_TEST(pMesh->getVertexDistribution().at(2).size() == 0);
      BOOST_TEST(pMesh->getVertexDistribution().at(0).at(0) == 0);
      BOOST_TEST(pMesh->getVertexDistribution().at(0).at(1) == 1);
      BOOST_TEST(pMesh->getVertexDistribution().at(1).at(0) == 1);
      BOOST_TEST(pMesh->vertices().size() == 2);
      BOOST_TEST(pMesh->vertices().at(0).getGlobalIndex() == 0);
      BOOST_TEST(pMesh->vertices().at(1).getGlobalIndex() == 1);
      BOOST_TEST(pMesh->vertices().at(0).isOwner() == true);
      BOOST_TEST(pMesh->vertices().at(1).isOwner() == false);
    } else if (context.isRank(1)) { //Slave2
      BOOST_TEST(pMesh->vertices().size() == 1);
      BOOST_TEST(pMesh->vertices().at(0).getGlobalIndex() == 1);
      BOOST_TEST(pMesh->vertices().at(0).isOwner() == true);
    } else if (context.isRank(2)) { //Slave3
      BOOST_TEST(pMesh->vertices().size() == 0);
    }
  }
}

BOOST_AUTO_TEST_CASE(ProvideAndReceiveCouplingMode)
{
  PRECICE_TEST("Fluid"_on(1_rank), "Solid"_on(1_rank), Require::Events);
  auto m2n = context.connectMasters("Solid", "Fluid");

  int  dimensions  = 2;
  bool flipNormals = false;

  if (context.isNamed("Solid")) {
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals, testing::nextMeshID()));
    createSolidzMesh2D(pSolidzMesh);
    ProvidedPartition part(pSolidzMesh);
    part.addM2N(m2n);
    part.communicate();
    part.compute();

    BOOST_TEST(pSolidzMesh->getGlobalNumberOfVertices() == 6);
    BOOST_TEST(pSolidzMesh->vertices().size() == 6);
    BOOST_TEST(pSolidzMesh->edges().size() == 5);
    BOOST_TEST(pSolidzMesh->getVertexOffsets().size() == 1);
    BOOST_TEST(pSolidzMesh->getVertexOffsets().at(0) == 6);
    BOOST_TEST(pSolidzMesh->vertices().at(0).getGlobalIndex() == 0);
    BOOST_TEST(pSolidzMesh->vertices().at(1).getGlobalIndex() == 1);
    BOOST_TEST(pSolidzMesh->vertices().at(2).getGlobalIndex() == 2);
    BOOST_TEST(pSolidzMesh->vertices().at(3).getGlobalIndex() == 3);
    BOOST_TEST(pSolidzMesh->vertices().at(4).getGlobalIndex() == 4);
    BOOST_TEST(pSolidzMesh->vertices().at(5).getGlobalIndex() == 5);
    BOOST_TEST(pSolidzMesh->vertices().at(0).isOwner() == true);
    BOOST_TEST(pSolidzMesh->vertices().at(1).isOwner() == true);
    BOOST_TEST(pSolidzMesh->vertices().at(2).isOwner() == true);
    BOOST_TEST(pSolidzMesh->vertices().at(3).isOwner() == true);
    BOOST_TEST(pSolidzMesh->vertices().at(4).isOwner() == true);
    BOOST_TEST(pSolidzMesh->vertices().at(5).isOwner() == true);
  } else {
    BOOST_TEST(context.isNamed("Fluid"));
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals, testing::nextMeshID()));

    mesh::PtrMesh       pOtherMesh(new mesh::Mesh("OtherMesh", dimensions, flipNormals, testing::nextMeshID()));
    mapping::PtrMapping boundingFromMapping = mapping::PtrMapping(
        new mapping::NearestNeighborMapping(mapping::Mapping::CONSISTENT, dimensions));
    boundingFromMapping->setMeshes(pSolidzMesh, pOtherMesh);

    double            safetyFactor = 0.1;
    ReceivedPartition part(pSolidzMesh, ReceivedPartition::ON_MASTER, safetyFactor);
    part.addFromMapping(boundingFromMapping);
    part.addM2N(m2n);
    part.communicate();
    part.compute();

    BOOST_TEST(pSolidzMesh->getGlobalNumberOfVertices() == 6);
    BOOST_TEST(pSolidzMesh->vertices().size() == 6);
    BOOST_TEST(pSolidzMesh->edges().size() == 5);
    BOOST_TEST(pSolidzMesh->getVertexOffsets().size() == 1);
    BOOST_TEST(pSolidzMesh->getVertexOffsets().at(0) == 6);
    BOOST_TEST(pSolidzMesh->vertices().at(0).getGlobalIndex() == 0);
    BOOST_TEST(pSolidzMesh->vertices().at(1).getGlobalIndex() == 1);
    BOOST_TEST(pSolidzMesh->vertices().at(2).getGlobalIndex() == 2);
    BOOST_TEST(pSolidzMesh->vertices().at(3).getGlobalIndex() == 3);
    BOOST_TEST(pSolidzMesh->vertices().at(4).getGlobalIndex() == 4);
    BOOST_TEST(pSolidzMesh->vertices().at(5).getGlobalIndex() == 5);
    BOOST_TEST(pSolidzMesh->vertices().at(0).isOwner() == true);
    BOOST_TEST(pSolidzMesh->vertices().at(1).isOwner() == true);
    BOOST_TEST(pSolidzMesh->vertices().at(2).isOwner() == true);
    BOOST_TEST(pSolidzMesh->vertices().at(3).isOwner() == true);
    BOOST_TEST(pSolidzMesh->vertices().at(4).isOwner() == true);
    BOOST_TEST(pSolidzMesh->vertices().at(5).isOwner() == true);
  }
}

BOOST_AUTO_TEST_CASE(TestCompareBoundingBoxes2D)
{
  PRECICE_TEST("SOLIDZ"_on(1_rank), "NASTIN"_on(3_ranks).setupMasterSlaves(), Require::Events);

  testing::ConnectionOptions options;
  options.useOnlyMasterCom = false;
  options.useTwoLevelInit  = true;
  auto m2n                 = context.connectMasters("SOLIDZ", "NASTIN", options);

  int  dimensions  = 2;
  bool flipNormals = true;

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
    std::vector<int>                connectedRanksList;
    int                             connectionMapSize = 0;
    std::map<int, std::vector<int>> receivedConnectionMap;
    mesh::PtrMesh                   pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals, testing::nextMeshID()));
    m2n->getMasterCommunication()->send(3, 0);
    com::CommunicateBoundingBox(m2n->getMasterCommunication()).sendBoundingBoxMap(sendGlobalBB, 0);
    m2n->getMasterCommunication()->receive(connectedRanksList, 0);
    connectionMapSize = connectedRanksList.size();
    BOOST_TEST_REQUIRE(connectionMapSize == 2);

    std::vector<int> connectedRanks;
    connectedRanks.push_back(-1);
    for (auto &rank : connectedRanksList) {
      receivedConnectionMap[rank] = connectedRanks;
    }

    com::CommunicateBoundingBox(m2n->getMasterCommunication()).receiveConnectionMap(receivedConnectionMap, 0);

    // test whether we receive correct connection map
    BOOST_TEST(receivedConnectionMap.at(0).at(0) == 2);
    BOOST_TEST(receivedConnectionMap.at(2).at(0) == 0);

  } else {
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals, testing::nextMeshID()));
    mesh::PtrMesh pNastinMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals, testing::nextMeshID()));

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
  tearDownParallelEnvironment();
}

BOOST_AUTO_TEST_CASE(TestCompareBoundingBoxes3D)
{
  PRECICE_TEST("SOLIDZ"_on(1_rank), "NASTIN"_on(3_ranks).setupMasterSlaves(), Require::Events);

  testing::ConnectionOptions options;
  options.useOnlyMasterCom = false;
  options.useTwoLevelInit  = true;
  auto m2n                 = context.connectMasters("SOLIDZ", "NASTIN", options);

  int  dimensions  = 3;
  bool flipNormals = true;

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
    std::vector<int>                connectedRanksList;
    int                             connectionMapSize = 0;
    std::map<int, std::vector<int>> receivedConnectionMap;
    mesh::PtrMesh                   pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals, testing::nextMeshID()));
    m2n->getMasterCommunication()->send(3, 0);
    com::CommunicateBoundingBox(m2n->getMasterCommunication()).sendBoundingBoxMap(sendGlobalBB, 0);
    m2n->getMasterCommunication()->receive(connectedRanksList, 0);
    connectionMapSize = connectedRanksList.size();
    BOOST_TEST(connectionMapSize == 2);

    std::vector<int> connectedRanks;
    connectedRanks.push_back(-1);
    for (auto &rank : connectedRanksList) {
      receivedConnectionMap[rank] = connectedRanks;
    }

    com::CommunicateBoundingBox(m2n->getMasterCommunication()).receiveConnectionMap(receivedConnectionMap, 0);

    // test whether we receive correct connection map
    BOOST_TEST(receivedConnectionMap.at(0).at(0) == 2);
    BOOST_TEST(receivedConnectionMap.at(2).at(0) == 0);

  } else {
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals, testing::nextMeshID()));
    mesh::PtrMesh pNastinMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals, testing::nextMeshID()));

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
  tearDownParallelEnvironment();
}

// Test with two "from" and two "to" mappings
BOOST_AUTO_TEST_CASE(RePartitionMultipleMappings)
{
  PRECICE_TEST("Solid"_on(1_rank), "Fluid"_on(3_ranks).setupMasterSlaves(), Require::Events);
  auto m2n = context.connectMasters("Solid", "Fluid");

  int             dimensions  = 2;
  bool            flipNormals = false;
  Eigen::VectorXd offset      = Eigen::VectorXd::Zero(dimensions);

  if (context.isNamed("Solid")) { //SOLIDZ
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals, testing::nextMeshID()));
    createSolidzMesh2D(pSolidzMesh);
    BOOST_TEST(pSolidzMesh->vertices().size() == 6);
    ProvidedPartition part(pSolidzMesh);
    part.addM2N(m2n);
    part.communicate();
  } else {
    BOOST_TEST(context.isNamed("Fluid"));
    mesh::PtrMesh pNastinMesh1(new mesh::Mesh("NastinMesh1", dimensions, flipNormals, testing::nextMeshID()));
    mesh::PtrMesh pNastinMesh2(new mesh::Mesh("NastinMesh1", dimensions, flipNormals, testing::nextMeshID()));
    mesh::PtrMesh pNastinMesh3(new mesh::Mesh("NastinMesh1", dimensions, flipNormals, testing::nextMeshID()));
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals, testing::nextMeshID()));

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
    pNastinMesh1->computeState();
    pNastinMesh1->computeBoundingBox();
    pNastinMesh2->computeState();
    pNastinMesh2->computeBoundingBox();
    pNastinMesh3->computeState();
    pNastinMesh3->computeBoundingBox();

    double safetyFactor = 0.1;

    ReceivedPartition part(pSolidzMesh, ReceivedPartition::ON_SLAVES, safetyFactor);
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
      if (context.isMaster()) { //Master
        BOOST_TEST(pSolidzMesh->vertices().size() == 2);
        BOOST_TEST(pSolidzMesh->edges().size() == 1);
      } else if (context.isRank(1)) { //Slave1
        BOOST_TEST(pSolidzMesh->vertices().size() == 0);
        BOOST_TEST(pSolidzMesh->edges().size() == 0);
      } else if (context.isRank(2)) { //Slave2
        BOOST_TEST(pSolidzMesh->vertices().size() == 2);
        BOOST_TEST(pSolidzMesh->edges().size() == 1);
      }
    }
  }

  tearDownParallelEnvironment();
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

#endif // PRECICE_NO_MPI
