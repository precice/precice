#include <Eigen/Core>
#include <algorithm>
#include <memory>
#include <ostream>
#include <string>
#include <utility>
#include <vector>
#include "logging/Logger.hpp"
#include "mapping/Mapping.hpp"
#include "mapping/RadialBasisFctMapping.hpp"
#include "mapping/config/MappingConfiguration.hpp"
#include "mapping/impl/BasisFunctions.hpp"
#include "mapping/tests/RadialBasisFctHelper.hpp"
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Utils.hpp"
#include "mesh/Vertex.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"

using namespace precice;
using namespace precice::mesh;
using namespace precice::mapping;
using namespace precice::testing;
using precice::testing::TestContext;

BOOST_AUTO_TEST_SUITE(MappingTests)
BOOST_AUTO_TEST_SUITE(RadialBasisFunctionMapping)

BOOST_AUTO_TEST_SUITE(Parallel)

/// Holds rank, owner, position and value of a single vertex
struct VertexSpecification {
  int                 rank;
  int                 owner;
  std::vector<double> position;
  std::vector<double> value;
};

/*
MeshSpecification format:
{ {rank, owner rank, {x, y, z}, {v}}, ... }

also see struct VertexSpecification.

- -1 on rank means all ranks
- -1 on owner rank means no rank
- x, y, z is position of vertex, z is optional, 2D mesh will be created then
- v is the value of the respective vertex. Only 1D supported at this time.

ReferenceSpecification format:
{ {rank, {v}, ... }
- -1 on rank means all ranks
- v is the expected value of n-th vertex on that particular rank
*/
using MeshSpecification = std::vector<VertexSpecification>;

/// Contains which values are expected on which rank: rank -> vector of data.
using ReferenceSpecification = std::vector<std::pair<int, std::vector<double>>>;

void getDistributedMesh(const TestContext &      context,
                        MeshSpecification const &vertices,
                        mesh::PtrMesh &          mesh,
                        mesh::PtrData &          data,
                        int                      globalIndexOffset = 0,
                        bool                     meshIsSmaller     = false)
{
  Eigen::VectorXd d;

  int i = 0;
  for (auto &vertex : vertices) {
    if (vertex.rank == context.rank or vertex.rank == -1) {
      if (vertex.position.size() == 3)      // 3-dimensional
        mesh->createVertex(Eigen::Vector3d(vertex.position.data()));
      else if (vertex.position.size() == 2) // 2-dimensional
        mesh->createVertex(Eigen::Vector2d(vertex.position.data()));

      int valueDimension = vertex.value.size();

      if (vertex.owner == context.rank)
        mesh->vertices().back().setOwner(true);
      else
        mesh->vertices().back().setOwner(false);

      d.conservativeResize(i * valueDimension + valueDimension);
      // Get data in every dimension
      for (int dim = 0; dim < valueDimension; ++dim) {
        d(i * valueDimension + dim) = vertex.value.at(dim);
      }
      i++;
    }
  }
  addGlobalIndex(mesh, globalIndexOffset);
  mesh->allocateDataValues();
  // All tests use eight vertices
  if (meshIsSmaller) {
    mesh->setGlobalNumberOfVertices(7);
  } else {
    mesh->setGlobalNumberOfVertices(8);
  }
  data->values() = d;
}

void testDistributed(const TestContext     &context,
                     Mapping               &mapping,
                     MeshSpecification      inMeshSpec,
                     MeshSpecification      outMeshSpec,
                     ReferenceSpecification referenceSpec,
                     int                    inGlobalIndexOffset = 0,
                     bool                   meshIsSmaller       = false)
{
  int meshDimension  = inMeshSpec.at(0).position.size();
  int valueDimension = inMeshSpec.at(0).value.size();

  mesh::PtrMesh inMesh(new mesh::Mesh("InMesh", meshDimension, testing::nextMeshID()));
  mesh::PtrData inData   = inMesh->createData("InData", valueDimension, 0_dataID);
  int           inDataID = inData->getID();

  getDistributedMesh(context, inMeshSpec, inMesh, inData, inGlobalIndexOffset);

  mesh::PtrMesh outMesh(new mesh::Mesh("outMesh", meshDimension, testing::nextMeshID()));
  mesh::PtrData outData   = outMesh->createData("OutData", valueDimension, 1_dataID);
  int           outDataID = outData->getID();

  getDistributedMesh(context, outMeshSpec, outMesh, outData, 0, meshIsSmaller);

  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  mapping.computeMapping();
  BOOST_TEST(mapping.hasComputedMapping() == true);
  mapping.map(inDataID, outDataID);

  int index = 0;
  for (auto &referenceVertex : referenceSpec) {
    if (referenceVertex.first == context.rank or referenceVertex.first == -1) {
      for (int dim = 0; dim < valueDimension; ++dim) {
        BOOST_TEST_INFO("Index of vertex: " << index << " - Dimension: " << dim);
        BOOST_TEST(outData->values()(index * valueDimension + dim) == referenceVertex.second.at(dim));
      }
      ++index;
    }
  }
  BOOST_TEST(outData->values().size() == index * valueDimension);
}

/// Test with a homogeneous distribution of mesh among ranks
BOOST_AUTO_TEST_CASE(DistributedConsistent2DV1)
{
  PRECICE_TEST(""_on(4_ranks).setupIntraComm());
  Multiquadrics fct(5.0);

  MeshSpecification in{// Consistent mapping: The inMesh is communicated
                       {-1, 0, {0, 0}, {1}},
                       {-1, 0, {0, 1}, {2}},
                       {-1, 1, {1, 0}, {3}},
                       {-1, 1, {1, 1}, {4}},
                       {-1, 2, {2, 0}, {5}},
                       {-1, 2, {2, 1}, {6}},
                       {-1, 3, {3, 0}, {7}},
                       {-1, 3, {3, 1}, {8}}};
  MeshSpecification out{// The outMesh is local, distributed among all ranks
                        {0, -1, {0, 0}, {0}},
                        {0, -1, {0, 1}, {0}},
                        {1, -1, {1, 0}, {0}},
                        {1, -1, {1, 1}, {0}},
                        {2, -1, {2, 0}, {0}},
                        {2, -1, {2, 1}, {0}},
                        {3, -1, {3, 0}, {0}},
                        {3, -1, {3, 1}, {0}}};

  ReferenceSpecification ref{// Tests for {0, 1} on the first rank, {1, 2} on the second, ...
                             {0, {1}},
                             {0, {2}},
                             {1, {3}},
                             {1, {4}},
                             {2, {5}},
                             {2, {6}},
                             {3, {7}},
                             {3, {8}}};

  RadialBasisFctMapping<RadialBasisFctSolver<Multiquadrics>> mapping_on(Mapping::CONSISTENT, 2, fct, {{false, false, false}}, Polynomial::ON);
  testDistributed(context, mapping_on, in, out, ref);
  RadialBasisFctMapping<RadialBasisFctSolver<Multiquadrics>> mapping_sep(Mapping::CONSISTENT, 2, fct, {{false, false, false}}, Polynomial::SEPARATE);
  testDistributed(context, mapping_sep, in, out, ref);
  RadialBasisFctMapping<RadialBasisFctSolver<Multiquadrics>> mapping_off(Mapping::CONSISTENT, 2, fct, {{false, false, false}}, Polynomial::OFF);
  testDistributed(context, mapping_off, in, out, ref);
}

BOOST_AUTO_TEST_CASE(DistributedConsistent2DV1Vector)
{
  PRECICE_TEST(""_on(4_ranks).setupIntraComm());
  Multiquadrics fct(5.0);

  MeshSpecification      in{// Consistent mapping: The inMesh is communicated
                       {-1, 0, {0, 0}, {1, 4}},
                       {-1, 0, {0, 1}, {2, 5}},
                       {-1, 1, {1, 0}, {3, 6}},
                       {-1, 1, {1, 1}, {4, 7}},
                       {-1, 2, {2, 0}, {5, 8}},
                       {-1, 2, {2, 1}, {6, 9}},
                       {-1, 3, {3, 0}, {7, 10}},
                       {-1, 3, {3, 1}, {8, 11}}};
  MeshSpecification      out{// The outMesh is local, distributed among all ranks
                        {0, -1, {0, 0}, {0, 0}},
                        {0, -1, {0, 1}, {0, 0}},
                        {1, -1, {1, 0}, {0, 0}},
                        {1, -1, {1, 1}, {0, 0}},
                        {2, -1, {2, 0}, {0, 0}},
                        {2, -1, {2, 1}, {0, 0}},
                        {3, -1, {3, 0}, {0, 0}},
                        {3, -1, {3, 1}, {0, 0}}};
  ReferenceSpecification ref{// Tests for {0, 1} on the first rank, {1, 2} on the second, ...
                             {0, {1, 4}},
                             {0, {2, 5}},
                             {1, {3, 6}},
                             {1, {4, 7}},
                             {2, {5, 8}},
                             {2, {6, 9}},
                             {3, {7, 10}},
                             {3, {8, 11}}};

  RadialBasisFctMapping<RadialBasisFctSolver<Multiquadrics>> mapping_on(Mapping::CONSISTENT, 2, fct, {{false, false, false}}, Polynomial::ON);
  testDistributed(context, mapping_on, in, out, ref);
  RadialBasisFctMapping<RadialBasisFctSolver<Multiquadrics>> mapping_sep(Mapping::CONSISTENT, 2, fct, {{false, false, false}}, Polynomial::SEPARATE);
  testDistributed(context, mapping_sep, in, out, ref);
  RadialBasisFctMapping<RadialBasisFctSolver<Multiquadrics>> mapping_off(Mapping::CONSISTENT, 2, fct, {{false, false, false}}, Polynomial::OFF);
  testDistributed(context, mapping_off, in, out, ref);
}

/// Using a more heterogeneous distributon of vertices and owner
BOOST_AUTO_TEST_CASE(DistributedConsistent2DV2)
{
  PRECICE_TEST(""_on(4_ranks).setupIntraComm());
  Multiquadrics fct(5.0);

  MeshSpecification                                          in{// Consistent mapping: The inMesh is communicated, rank 2 owns no vertices
                       {-1, 0, {0, 0}, {1}},
                       {-1, 0, {0, 1}, {2}},
                       {-1, 1, {1, 0}, {3}},
                       {-1, 1, {1, 1}, {4}},
                       {-1, 1, {2, 0}, {5}},
                       {-1, 3, {2, 1}, {6}},
                       {-1, 3, {3, 0}, {7}},
                       {-1, 3, {3, 1}, {8}}};
  MeshSpecification                                          out{// The outMesh is local, rank 1 is empty
                        {0, -1, {0, 0}, {0}},
                        {0, -1, {0, 1}, {0}},
                        {0, -1, {1, 0}, {0}},
                        {2, -1, {1, 1}, {0}},
                        {2, -1, {2, 0}, {0}},
                        {2, -1, {2, 1}, {0}},
                        {3, -1, {3, 0}, {0}},
                        {3, -1, {3, 1}, {0}}};
  ReferenceSpecification                                     ref{// Tests for {0, 1, 2} on the first rank,
                             // second rank (consistent with the outMesh) is empty, ...
                             {0, {1}},
                             {0, {2}},
                             {0, {3}},
                             {2, {4}},
                             {2, {5}},
                             {2, {6}},
                             {3, {7}},
                             {3, {8}}};
  RadialBasisFctMapping<RadialBasisFctSolver<Multiquadrics>> mapping_on(Mapping::CONSISTENT, 2, fct, {{false, false, false}}, Polynomial::ON);
  testDistributed(context, mapping_on, in, out, ref);
  RadialBasisFctMapping<RadialBasisFctSolver<Multiquadrics>> mapping_sep(Mapping::CONSISTENT, 2, fct, {{false, false, false}}, Polynomial::SEPARATE);
  testDistributed(context, mapping_sep, in, out, ref);
  RadialBasisFctMapping<RadialBasisFctSolver<Multiquadrics>> mapping_off(Mapping::CONSISTENT, 2, fct, {{false, false, false}}, Polynomial::OFF);
  testDistributed(context, mapping_off, in, out, ref);
}

/// Test with a very heterogeneous distributed and non-continuous ownership
BOOST_AUTO_TEST_CASE(DistributedConsistent2DV3)
{
  PRECICE_TEST(""_on(4_ranks).setupIntraComm());
  Multiquadrics fct(5.0);

  std::vector<int> globalIndexOffsets = {0, 0, 0, 4};

  MeshSpecification in{
      // Rank 0 has part of the mesh, owns a subpart
      {0, 0, {0, 0}, {1}},
      {0, 0, {0, 1}, {2}},
      {0, 0, {1, 0}, {3}},
      {0, -1, {1, 1}, {4}},
      {0, -1, {2, 0}, {5}},
      {0, -1, {2, 1}, {6}},
      // Rank 1 has no vertices
      // Rank 2 has the entire mesh, but owns just 3 and 5.
      {2, -1, {0, 0}, {1}},
      {2, -1, {0, 1}, {2}},
      {2, -1, {1, 0}, {3}},
      {2, 2, {1, 1}, {4}},
      {2, -1, {2, 0}, {5}},
      {2, 2, {2, 1}, {6}},
      {2, -1, {3, 0}, {7}},
      {2, -1, {3, 1}, {8}},
      // Rank 3 has the last 4 vertices, owns 4, 6 and 7
      {3, 3, {2, 0}, {5}},
      {3, -1, {2, 1}, {6}},
      {3, 3, {3, 0}, {7}},
      {3, 3, {3, 1}, {8}},
  };
  MeshSpecification      out{// The outMesh is local, rank 1 is empty
                        {0, -1, {0, 0}, {0}},
                        {0, -1, {0, 1}, {0}},
                        {0, -1, {1, 0}, {0}},
                        {2, -1, {1, 1}, {0}},
                        {2, -1, {2, 0}, {0}},
                        {2, -1, {2, 1}, {0}},
                        {3, -1, {3, 0}, {0}},
                        {3, -1, {3, 1}, {0}}};
  ReferenceSpecification ref{// Tests for {0, 1, 2} on the first rank,
                             // second rank (consistent with the outMesh) is empty, ...
                             {0, {1}},
                             {0, {2}},
                             {0, {3}},
                             {2, {4}},
                             {2, {5}},
                             {2, {6}},
                             {3, {7}},
                             {3, {8}}};

  RadialBasisFctMapping<RadialBasisFctSolver<Multiquadrics>> mapping_on(Mapping::CONSISTENT, 2, fct, {{false, false, false}}, Polynomial::ON);
  testDistributed(context, mapping_on, in, out, ref, globalIndexOffsets.at(context.rank));
  RadialBasisFctMapping<RadialBasisFctSolver<Multiquadrics>> mapping_sep(Mapping::CONSISTENT, 2, fct, {{false, false, false}}, Polynomial::SEPARATE);
  testDistributed(context, mapping_sep, in, out, ref, globalIndexOffsets.at(context.rank));
  RadialBasisFctMapping<RadialBasisFctSolver<Multiquadrics>> mapping_off(Mapping::CONSISTENT, 2, fct, {{false, false, false}}, Polynomial::OFF);
  testDistributed(context, mapping_off, in, out, ref, globalIndexOffsets.at(context.rank));
}

/// Test with a very heterogeneous distributed and non-continuous ownership
BOOST_AUTO_TEST_CASE(DistributedConsistent2DV3Vector)
{
  PRECICE_TEST(""_on(4_ranks).setupIntraComm());
  Multiquadrics fct(5.0);

  std::vector<int> globalIndexOffsets = {0, 0, 0, 4};

  MeshSpecification in{
      // Rank 0 has part of the mesh, owns a subpart
      {0, 0, {0, 0}, {1, 4}},
      {0, 0, {0, 1}, {2, 5}},
      {0, 0, {1, 0}, {3, 6}},
      {0, -1, {1, 1}, {4, 7}},
      {0, -1, {2, 0}, {5, 8}},
      {0, -1, {2, 1}, {6, 9}},
      // Rank 1 has no vertices
      // Rank 2 has the entire mesh, but owns just 3 and 5.
      {2, -1, {0, 0}, {1, 4}},
      {2, -1, {0, 1}, {2, 5}},
      {2, -1, {1, 0}, {3, 6}},
      {2, 2, {1, 1}, {4, 7}},
      {2, -1, {2, 0}, {5, 8}},
      {2, 2, {2, 1}, {6, 9}},
      {2, -1, {3, 0}, {7, 10}},
      {2, -1, {3, 1}, {8, 11}},
      // Rank 3 has the last 4 vertices, owns 4, 6 and 7
      {3, 3, {2, 0}, {5, 8}},
      {3, -1, {2, 1}, {6, 9}},
      {3, 3, {3, 0}, {7, 10}},
      {3, 3, {3, 1}, {8, 11}},
  };
  MeshSpecification                                          out{// The outMesh is local, rank 1 is empty
                        {0, -1, {0, 0}, {0, 0}},
                        {0, -1, {0, 1}, {0, 0}},
                        {0, -1, {1, 0}, {0, 0}},
                        {2, -1, {1, 1}, {0, 0}},
                        {2, -1, {2, 0}, {0, 0}},
                        {2, -1, {2, 1}, {0, 0}},
                        {3, -1, {3, 0}, {0, 0}},
                        {3, -1, {3, 1}, {0, 0}}};
  ReferenceSpecification                                     ref{// Tests for {0, 1, 2} on the first rank,
                             // second rank (consistent with the outMesh) is empty, ...
                             {0, {1, 4}},
                             {0, {2, 5}},
                             {0, {3, 6}},
                             {2, {4, 7}},
                             {2, {5, 8}},
                             {2, {6, 9}},
                             {3, {7, 10}},
                             {3, {8, 11}}};
  RadialBasisFctMapping<RadialBasisFctSolver<Multiquadrics>> mapping_on(Mapping::CONSISTENT, 2, fct, {{false, false, false}}, Polynomial::ON);
  testDistributed(context, mapping_on, in, out, ref, globalIndexOffsets.at(context.rank));
  RadialBasisFctMapping<RadialBasisFctSolver<Multiquadrics>> mapping_sep(Mapping::CONSISTENT, 2, fct, {{false, false, false}}, Polynomial::SEPARATE);
  testDistributed(context, mapping_sep, in, out, ref, globalIndexOffsets.at(context.rank));
  RadialBasisFctMapping<RadialBasisFctSolver<Multiquadrics>> mapping_off(Mapping::CONSISTENT, 2, fct, {{false, false, false}}, Polynomial::OFF);
  testDistributed(context, mapping_off, in, out, ref, globalIndexOffsets.at(context.rank));
}

/// Some ranks are empty, does not converge
BOOST_AUTO_TEST_CASE(DistributedConsistent2DV4)
{
  PRECICE_TEST(""_on(4_ranks).setupIntraComm());
  ThinPlateSplines fct;

  std::vector<int> globalIndexOffsets = {0, 0, 0, 0};

  MeshSpecification in{
      // Rank 0 has no vertices
      // Rank 1 has the entire mesh, owns a subpart
      {1, 1, {0, 0}, {1.1}},
      {1, 1, {0, 1}, {2.5}},
      {1, 1, {1, 0}, {3}},
      {1, 1, {1, 1}, {4}},
      {1, -1, {2, 0}, {5}},
      {1, -1, {2, 1}, {6}},
      {1, -1, {3, 0}, {7}},
      {1, -1, {3, 1}, {8}},
      // Rank 2 has the entire mesh, owns a subpart
      {2, -1, {0, 0}, {1.1}},
      {2, -1, {0, 1}, {2.5}},
      {2, -1, {1, 0}, {3}},
      {2, -1, {1, 1}, {4}},
      {2, 2, {2, 0}, {5}},
      {2, 2, {2, 1}, {6}},
      {2, 2, {3, 0}, {7}},
      {2, 2, {3, 1}, {8}},
      // Rank 3 has no vertices
  };
  MeshSpecification                                             out{// The outMesh is local, rank 0 and 3 are empty
                        // not same order as input mesh and vertex (2,0) appears twice
                        {1, -1, {2, 0}, {0}},
                        {1, -1, {1, 0}, {0}},
                        {1, -1, {0, 1}, {0}},
                        {1, -1, {1, 1}, {0}},
                        {1, -1, {0, 0}, {0}},
                        {2, -1, {2, 0}, {0}},
                        {2, -1, {2, 1}, {0}},
                        {2, -1, {3, 0}, {0}},
                        {2, -1, {3, 1}, {0}}};
  ReferenceSpecification ref{{1, {5}},
                                                                         {1, {3}},
                                                                         {1, {2.5}},
                                                                         {1, {4}},
                                                                         {1, {1.1}},
                                                                         {2, {5}},
                                                                         {2, {6}},
                                                                         {2, {7}},
                                                                         {2, {8}}};
  RadialBasisFctMapping<RadialBasisFctSolver<ThinPlateSplines>> mapping_on(Mapping::CONSISTENT, 2, fct, {{false, false, false}}, Polynomial::ON);
  testDistributed(context, mapping_on, in, out, ref, globalIndexOffsets.at(context.rank));
  RadialBasisFctMapping<RadialBasisFctSolver<ThinPlateSplines>> mapping_sep(Mapping::CONSISTENT, 2, fct, {{false, false, false}}, Polynomial::SEPARATE);
  testDistributed(context, mapping_sep, in, out, ref, globalIndexOffsets.at(context.rank));
  RadialBasisFctMapping<RadialBasisFctSolver<ThinPlateSplines>> mapping_off(Mapping::CONSISTENT, 2, fct, {{false, false, false}}, Polynomial::OFF);
  testDistributed(context, mapping_off, in, out, ref, globalIndexOffsets.at(context.rank));
}

// same as 2DV4, but all ranks have vertices
BOOST_AUTO_TEST_CASE(DistributedConsistent2DV5)
{
  PRECICE_TEST(""_on(4_ranks).setupIntraComm());
  ThinPlateSplines fct;

  std::vector<int>  globalIndexOffsets = {0, 0, 0, 0};
  MeshSpecification in{
      // Every rank has the entire mesh and owns a subpart
      {0, 0, {0, 0}, {1.1}},
      {0, 0, {0, 1}, {2.5}},
      {0, -1, {1, 0}, {3}},
      {0, -1, {1, 1}, {4}},
      {0, -1, {2, 0}, {5}},
      {0, -1, {2, 1}, {6}},
      {0, -1, {3, 0}, {7}},
      {0, -1, {3, 1}, {8}},
      {1, -1, {0, 0}, {1.1}},
      {1, -1, {0, 1}, {2.5}},
      {1, 1, {1, 0}, {3}},
      {1, 1, {1, 1}, {4}},
      {1, -1, {2, 0}, {5}},
      {1, -1, {2, 1}, {6}},
      {1, -1, {3, 0}, {7}},
      {1, -1, {3, 1}, {8}},
      {2, -1, {0, 0}, {1.1}},
      {2, -1, {0, 1}, {2.5}},
      {2, -1, {1, 0}, {3}},
      {2, -1, {1, 1}, {4}},
      {2, 2, {2, 0}, {5}},
      {2, 2, {2, 1}, {6}},
      {2, -1, {3, 0}, {7}},
      {2, -1, {3, 1}, {8}},
      {3, -1, {0, 0}, {1.1}},
      {3, -1, {0, 1}, {2.5}},
      {3, -1, {1, 0}, {3}},
      {3, -1, {1, 1}, {4}},
      {3, -1, {2, 0}, {5}},
      {3, -1, {2, 1}, {6}},
      {3, 3, {3, 0}, {7}},
      {3, 3, {3, 1}, {8}},
  };
  MeshSpecification                                             out{// The outMesh is local, rank 0 and 3 are empty
                        // not same order as input mesh and vertex (2,0) appears twice
                        {1, -1, {2, 0}, {0}},
                        {1, -1, {1, 0}, {0}},
                        {1, -1, {0, 1}, {0}},
                        {1, -1, {1, 1}, {0}},
                        {1, -1, {0, 0}, {0}},
                        {2, -1, {2, 0}, {0}},
                        {2, -1, {2, 1}, {0}},
                        {2, -1, {3, 0}, {0}},
                        {2, -1, {3, 1}, {0}}};
  ReferenceSpecification ref{{1, {5}},
                                                                         {1, {3}},
                                                                         {1, {2.5}},
                                                                         {1, {4}},
                                                                         {1, {1.1}},
                                                                         {2, {5}},
                                                                         {2, {6}},
                                                                         {2, {7}},
                                                                         {2, {8}}};
  RadialBasisFctMapping<RadialBasisFctSolver<ThinPlateSplines>> mapping_on(Mapping::CONSISTENT, 2, fct, {{false, false, false}}, Polynomial::ON);
  testDistributed(context, mapping_on, in, out, ref, globalIndexOffsets.at(context.rank));
  RadialBasisFctMapping<RadialBasisFctSolver<ThinPlateSplines>> mapping_sep(Mapping::CONSISTENT, 2, fct, {{false, false, false}}, Polynomial::SEPARATE);
  testDistributed(context, mapping_sep, in, out, ref, globalIndexOffsets.at(context.rank));
  RadialBasisFctMapping<RadialBasisFctSolver<ThinPlateSplines>> mapping_off(Mapping::CONSISTENT, 2, fct, {{false, false, false}}, Polynomial::OFF);
  testDistributed(context, mapping_off, in, out, ref, globalIndexOffsets.at(context.rank));
}

/// same as 2DV4, but strictly linear input values, converges and gives correct results
BOOST_AUTO_TEST_CASE(DistributedConsistent2DV6,
                     *boost::unit_test::tolerance(1e-7))
{
  PRECICE_TEST(""_on(4_ranks).setupIntraComm());
  ThinPlateSplines fct;
  std::vector<int> globalIndexOffsets = {0, 0, 0, 0};

  MeshSpecification in{
      // Rank 0 has no vertices
      // Rank 1 has the entire mesh, owns a subpart
      {1, 1, {0, 0}, {1}},
      {1, 1, {0, 1}, {2}},
      {1, 1, {1, 0}, {3}},
      {1, 1, {1, 1}, {4}},
      {1, -1, {2, 0}, {5}},
      {1, -1, {2, 1}, {6}},
      {1, -1, {3, 0}, {7}},
      {1, -1, {3, 1}, {8}},
      // Rank 2 has the entire mesh, owns a subpart
      {2, -1, {0, 0}, {1}},
      {2, -1, {0, 1}, {2}},
      {2, -1, {1, 0}, {3}},
      {2, -1, {1, 1}, {4}},
      {2, 2, {2, 0}, {5}},
      {2, 2, {2, 1}, {6}},
      {2, 2, {3, 0}, {7}},
      {2, 2, {3, 1}, {8}},
      // Rank 3 has no vertices
  };
  MeshSpecification                                             out{// The outMesh is local, rank 0 and 3 are empty
                        // not same order as input mesh and vertex (2,0) appears twice
                        {1, -1, {2, 0}, {0}},
                        {1, -1, {1, 0}, {0}},
                        {1, -1, {0, 1}, {0}},
                        {1, -1, {1, 1}, {0}},
                        {1, -1, {0, 0}, {0}},
                        {2, -1, {2, 0}, {0}},
                        {2, -1, {2, 1}, {0}},
                        {2, -1, {3, 0}, {0}},
                        {2, -1, {3, 1}, {0}}};
  ReferenceSpecification ref{{1, {5}},
                                                                         {1, {3}},
                                                                         {1, {2}},
                                                                         {1, {4}},
                                                                         {1, {1}},
                                                                         {2, {5}},
                                                                         {2, {6}},
                                                                         {2, {7}},
                                                                         {2, {8}}};
  RadialBasisFctMapping<RadialBasisFctSolver<ThinPlateSplines>> mapping_on(Mapping::CONSISTENT, 2, fct, {{false, false, false}}, Polynomial::ON);
  testDistributed(context, mapping_on, in, out, ref, globalIndexOffsets.at(context.rank));
  RadialBasisFctMapping<RadialBasisFctSolver<ThinPlateSplines>> mapping_sep(Mapping::CONSISTENT, 2, fct, {{false, false, false}}, Polynomial::SEPARATE);
  testDistributed(context, mapping_sep, in, out, ref, globalIndexOffsets.at(context.rank));
  RadialBasisFctMapping<RadialBasisFctSolver<ThinPlateSplines>> mapping_off(Mapping::CONSISTENT, 2, fct, {{false, false, false}}, Polynomial::OFF);
  testDistributed(context, mapping_off, in, out, ref, globalIndexOffsets.at(context.rank));
}

/// Test with a homogeneous distribution of mesh among ranks
BOOST_AUTO_TEST_CASE(DistributedConservative2DV1)
{
  PRECICE_TEST(""_on(4_ranks).setupIntraComm());
  Multiquadrics fct(5.0);

  MeshSpecification      in{// Conservative mapping: The inMesh is local
                       {0, -1, {0, 0}, {1}},
                       {0, -1, {0, 1}, {2}},
                       {1, -1, {1, 0}, {3}},
                       {1, -1, {1, 1}, {4}},
                       {2, -1, {2, 0}, {5}},
                       {2, -1, {2, 1}, {6}},
                       {3, -1, {3, 0}, {7}},
                       {3, -1, {3, 1}, {8}}};
  MeshSpecification      out{// The outMesh is distributed
                        {-1, 0, {0, 0}, {0}},
                        {-1, 0, {0, 1}, {0}},
                        {-1, 1, {1, 0}, {0}},
                        {-1, 1, {1, 1}, {0}},
                        {-1, 2, {2, 0}, {0}},
                        {-1, 2, {2, 1}, {0}},
                        {-1, 3, {3, 0}, {0}},
                        {-1, 3, {3, 1}, {0}}};
  ReferenceSpecification ref{// Tests for {0, 1, 0, 0, 0, 0, 0, 0} on the first rank,
                             // {0, 0, 2, 3, 0, 0, 0, 0} on the second, ...
                             {0, {1}},
                             {0, {2}},
                             {0, {0}},
                             {0, {0}},
                             {0, {0}},
                             {0, {0}},
                             {0, {0}},
                             {0, {0}},
                             {1, {0}},
                             {1, {0}},
                             {1, {3}},
                             {1, {4}},
                             {1, {0}},
                             {1, {0}},
                             {1, {0}},
                             {1, {0}},
                             {2, {0}},
                             {2, {0}},
                             {2, {0}},
                             {2, {0}},
                             {2, {5}},
                             {2, {6}},
                             {2, {0}},
                             {2, {0}},
                             {3, {0}},
                             {3, {0}},
                             {3, {0}},
                             {3, {0}},
                             {3, {0}},
                             {3, {0}},
                             {3, {7}},
                             {3, {8}}};

  RadialBasisFctMapping<RadialBasisFctSolver<Multiquadrics>> mapping_on(Mapping::CONSERVATIVE, 2, fct, {{false, false, false}}, Polynomial::ON);
  testDistributed(context, mapping_on, in, out, ref, context.rank * 2);
  RadialBasisFctMapping<RadialBasisFctSolver<Multiquadrics>> mapping_sep(Mapping::CONSERVATIVE, 2, fct, {{false, false, false}}, Polynomial::SEPARATE);
  testDistributed(context, mapping_sep, in, out, ref, context.rank * 2);
  RadialBasisFctMapping<RadialBasisFctSolver<Multiquadrics>> mapping_off(Mapping::CONSERVATIVE, 2, fct, {{false, false, false}}, Polynomial::OFF);
  testDistributed(context, mapping_off, in, out, ref, context.rank * 2);
}

/// Test with a homogeneous distribution of mesh among ranks
BOOST_AUTO_TEST_CASE(DistributedConservative2DV1Vector)
{
  PRECICE_TEST(""_on(4_ranks).setupIntraComm());
  Multiquadrics                                              fct(5.0);
  MeshSpecification                                          in{// Conservative mapping: The inMesh is local
                       {0, -1, {0, 0}, {1, 4}},
                       {0, -1, {0, 1}, {2, 5}},
                       {1, -1, {1, 0}, {3, 6}},
                       {1, -1, {1, 1}, {4, 7}},
                       {2, -1, {2, 0}, {5, 8}},
                       {2, -1, {2, 1}, {6, 9}},
                       {3, -1, {3, 0}, {7, 10}},
                       {3, -1, {3, 1}, {8, 11}}};
  MeshSpecification                                          out{// The outMesh is distributed
                        {-1, 0, {0, 0}, {0, 0}},
                        {-1, 0, {0, 1}, {0, 0}},
                        {-1, 1, {1, 0}, {0, 0}},
                        {-1, 1, {1, 1}, {0, 0}},
                        {-1, 2, {2, 0}, {0, 0}},
                        {-1, 2, {2, 1}, {0, 0}},
                        {-1, 3, {3, 0}, {0, 0}},
                        {-1, 3, {3, 1}, {0, 0}}};
  ReferenceSpecification                                     ref{// Tests for {0, 1, 0, 0, 0, 0, 0, 0} on the first rank,
                             // {0, 0, 2, 3, 0, 0, 0, 0} on the second, ...
                             {0, {1, 4}},
                             {0, {2, 5}},
                             {0, {0, 0}},
                             {0, {0, 0}},
                             {0, {0, 0}},
                             {0, {0, 0}},
                             {0, {0, 0}},
                             {0, {0, 0}},
                             {1, {0, 0}},
                             {1, {0, 0}},
                             {1, {3, 6}},
                             {1, {4, 7}},
                             {1, {0, 0}},
                             {1, {0, 0}},
                             {1, {0, 0}},
                             {1, {0, 0}},
                             {2, {0, 0}},
                             {2, {0, 0}},
                             {2, {0, 0}},
                             {2, {0, 0}},
                             {2, {5, 8}},
                             {2, {6, 9}},
                             {2, {0, 0}},
                             {2, {0, 0}},
                             {3, {0, 0}},
                             {3, {0, 0}},
                             {3, {0, 0}},
                             {3, {0, 0}},
                             {3, {0, 0}},
                             {3, {0, 0}},
                             {3, {7, 10}},
                             {3, {8, 11}}};
  RadialBasisFctMapping<RadialBasisFctSolver<Multiquadrics>> mapping_on(Mapping::CONSERVATIVE, 2, fct, {{false, false, false}}, Polynomial::ON);
  testDistributed(context, mapping_on, in, out, ref, context.rank * 2);
  RadialBasisFctMapping<RadialBasisFctSolver<Multiquadrics>> mapping_sep(Mapping::CONSERVATIVE, 2, fct, {{false, false, false}}, Polynomial::SEPARATE);
  testDistributed(context, mapping_sep, in, out, ref, context.rank * 2);
  RadialBasisFctMapping<RadialBasisFctSolver<Multiquadrics>> mapping_off(Mapping::CONSERVATIVE, 2, fct, {{false, false, false}}, Polynomial::OFF);
  testDistributed(context, mapping_off, in, out, ref, context.rank * 2);
}

/// Using a more heterogeneous distribution of vertices and owner
BOOST_AUTO_TEST_CASE(DistributedConservative2DV2)
{
  PRECICE_TEST(""_on(4_ranks).setupIntraComm())
  Multiquadrics fct(5.0);

  std::vector<int> globalIndexOffsets = {0, 0, 4, 6};

  MeshSpecification      in{// Conservative mapping: The inMesh is local but rank 0 has no vertices
                       {1, -1, {0, 0}, {1}},
                       {1, -1, {0, 1}, {2}},
                       {1, -1, {1, 0}, {3}},
                       {1, -1, {1, 1}, {4}},
                       {2, -1, {2, 0}, {5}},
                       {2, -1, {2, 1}, {6}},
                       {3, -1, {3, 0}, {7}},
                       {3, -1, {3, 1}, {8}}};
  MeshSpecification      out{// The outMesh is distributed, rank 0 owns no vertex
                        {-1, 1, {0, 0}, {0}},
                        {-1, 1, {0, 1}, {0}},
                        {-1, 1, {1, 0}, {0}},
                        {-1, 1, {1, 1}, {0}},
                        {-1, 2, {2, 0}, {0}},
                        {-1, 2, {2, 1}, {0}},
                        {-1, 3, {3, 0}, {0}},
                        {-1, 3, {3, 1}, {0}}};
  ReferenceSpecification ref{// Tests for {0, 0, 0, 0, 0, 0, 0, 0} on the first rank,
                             // {1, 2, 2, 3, 0, 0, 0, 0} on the second, ...
                             {0, {0}},
                             {0, {0}},
                             {0, {0}},
                             {0, {0}},
                             {0, {0}},
                             {0, {0}},
                             {0, {0}},
                             {0, {0}},
                             {1, {1}},
                             {1, {2}},
                             {1, {3}},
                             {1, {4}},
                             {1, {0}},
                             {1, {0}},
                             {1, {0}},
                             {1, {0}},
                             {2, {0}},
                             {2, {0}},
                             {2, {0}},
                             {2, {0}},
                             {2, {5}},
                             {2, {6}},
                             {2, {0}},
                             {2, {0}},
                             {3, {0}},
                             {3, {0}},
                             {3, {0}},
                             {3, {0}},
                             {3, {0}},
                             {3, {0}},
                             {3, {7}},
                             {3, {8}}};

  RadialBasisFctMapping<RadialBasisFctSolver<Multiquadrics>> mapping_on(Mapping::CONSERVATIVE, 2, fct, {{false, false, false}}, Polynomial::ON);
  testDistributed(context, mapping_on, in, out, ref, globalIndexOffsets.at(context.rank));
  RadialBasisFctMapping<RadialBasisFctSolver<Multiquadrics>> mapping_sep(Mapping::CONSERVATIVE, 2, fct, {{false, false, false}}, Polynomial::SEPARATE);
  testDistributed(context, mapping_sep, in, out, ref, globalIndexOffsets.at(context.rank));
  RadialBasisFctMapping<RadialBasisFctSolver<Multiquadrics>> mapping_off(Mapping::CONSERVATIVE, 2, fct, {{false, false, false}}, Polynomial::OFF);
  testDistributed(context, mapping_off, in, out, ref, globalIndexOffsets.at(context.rank));
}

/// Using meshes of different sizes, inMesh is smaller then outMesh
BOOST_AUTO_TEST_CASE(DistributedConservative2DV3)
{
  PRECICE_TEST(""_on(4_ranks).setupIntraComm());
  Multiquadrics    fct(2.0);
  std::vector<int> globalIndexOffsets = {0, 0, 3, 5};

  MeshSpecification      in{// Conservative mapping: The inMesh is local but rank 0 has no vertices
                       {1, -1, {0, 0}, {1}},
                       {1, -1, {1, 0}, {3}},
                       {1, -1, {1, 1}, {4}},
                       {2, -1, {2, 0}, {5}},
                       {2, -1, {2, 1}, {6}},
                       {3, -1, {3, 0}, {7}},
                       {3, -1, {3, 1}, {8}}}; // Sum of all vertices is 34
  MeshSpecification      out{                      // The outMesh is distributed, rank 0 owns no vertex
                        {-1, 1, {0, 0}, {0}},
                        {-1, 1, {0, 1}, {0}},
                        {-1, 1, {1, 0}, {0}},
                        {-1, 1, {1, 1}, {0}},
                        {-1, 2, {2, 0}, {0}},
                        {-1, 2, {2, 1}, {0}},
                        {-1, 3, {3, 0}, {0}},
                        {-1, 3, {3, 1}, {0}}};
  ReferenceSpecification ref{// Tests for {0, 0, 0, 0, 0, 0, 0, 0} on the first rank,
                             // {1, 2, 2, 3, 0, 0, 0, 0} on the second, ...
                             {0, {0}},
                             {0, {0}},
                             {0, {0}},
                             {0, {0}},
                             {0, {0}},
                             {0, {0}},
                             {0, {0}},
                             {0, {0}},
                             {1, {1}},
                             {1, {0}},
                             {1, {3}},
                             {1, {4}},
                             {1, {0}},
                             {1, {0}},
                             {1, {0}},
                             {1, {0}},
                             {2, {0}},
                             {2, {0}},
                             {2, {0}},
                             {2, {0}},
                             {2, {5}},
                             {2, {6}},
                             {2, {0}},
                             {2, {0}},
                             {3, {0}},
                             {3, {0}},
                             {3, {0}},
                             {3, {0}},
                             {3, {0}},
                             {3, {0}},
                             {3, {7}},
                             {3, {8}}};
  // Sum of reference is also 34
  RadialBasisFctMapping<RadialBasisFctSolver<Multiquadrics>> mapping_on(Mapping::CONSERVATIVE, 2, fct, {{false, false, false}}, Polynomial::ON);
  testDistributed(context, mapping_on, in, out, ref, globalIndexOffsets.at(context.rank));
  RadialBasisFctMapping<RadialBasisFctSolver<Multiquadrics>> mapping_sep(Mapping::CONSERVATIVE, 2, fct, {{false, false, false}}, Polynomial::SEPARATE);
  testDistributed(context, mapping_sep, in, out, ref, globalIndexOffsets.at(context.rank));
  RadialBasisFctMapping<RadialBasisFctSolver<Multiquadrics>> mapping_off(Mapping::CONSERVATIVE, 2, fct, {{false, false, false}}, Polynomial::OFF);
  testDistributed(context, mapping_off, in, out, ref, globalIndexOffsets.at(context.rank));
}

/// Using meshes of different sizes, outMesh is smaller then inMesh
BOOST_AUTO_TEST_CASE(DistributedConservative2DV4,
                     *boost::unit_test::tolerance(1e-6))
{
  PRECICE_TEST(""_on(4_ranks).setupIntraComm());
  ThinPlateSplines fct;
  std::vector<int> globalIndexOffsets = {0, 2, 4, 6};

  MeshSpecification                                             in{// Conservative mapping: The inMesh is local
                       {0, -1, {0, 0}, {1}},
                       {0, -1, {0, 1}, {2}},
                       {1, -1, {1, 0}, {3}},
                       {1, -1, {1, 1}, {4}},
                       {2, -1, {2, 0}, {5}},
                       {2, -1, {2, 1}, {6}},
                       {3, -1, {3, 0}, {7}},
                       {3, -1, {3, 1}, {8}}}; // Sum is 36
  MeshSpecification                                             out{                      // The outMesh is distributed, rank 0 has no vertex at all
                        {-1, 1, {0, 1}, {0}},
                        {-1, 1, {1, 0}, {0}},
                        {-1, 1, {1, 1}, {0}},
                        {-1, 2, {2, 0}, {0}},
                        {-1, 2, {2, 1}, {0}},
                        {-1, 3, {3, 0}, {0}},
                        {-1, 3, {3, 1}, {0}}};
  ReferenceSpecification                                        ref{// Tests for {0, 0, 0, 0, 0, 0, 0, 0} on the first rank,
                             // {2, 3, 4, 3, 0, 0, 0, 0} on the second, ...
                             {0, {0}},
                             {0, {0}},
                             {0, {0}},
                             {0, {0}},
                             {0, {0}},
                             {0, {0}},
                             {0, {0}},
                             {1, {3.1307987056295525}},
                             {1, {4.0734031184906971}},
                             {1, {3.0620533966233006}},
                             {1, {0}},
                             {1, {0}},
                             {1, {0}},
                             {1, {0}},
                             {2, {0}},
                             {2, {0}},
                             {2, {0}},
                             {2, {4.4582564956060873}},
                             {2, {5.8784343572772633}},
                             {2, {0}},
                             {2, {0}},
                             {3, {0}},
                             {3, {0}},
                             {3, {0}},
                             {3, {0}},
                             {3, {0}},
                             {3, {7.4683403859032156}},
                             {3, {7.9287135404698859}}}; // Sum is ~36
  RadialBasisFctMapping<RadialBasisFctSolver<ThinPlateSplines>> mapping_sep(Mapping::CONSERVATIVE, 2, fct, {{false, false, false}}, Polynomial::SEPARATE);
  testDistributed(context, mapping_sep, in, out, ref, globalIndexOffsets.at(context.rank), true);
  // Polynomial == OFF won't reach the desired accuracy
}

/// Tests a non-contigous owner distributed at the outMesh
BOOST_AUTO_TEST_CASE(testDistributedConservative2DV5)
{
  PRECICE_TEST(""_on(4_ranks).setupIntraComm());
  Multiquadrics          fct(5.0);
  MeshSpecification      in{// Conservative mapping: The inMesh is local
                       {0, -1, {0, 0}, {1}},
                       {0, -1, {0, 1}, {2}},
                       {1, -1, {1, 0}, {3}},
                       {1, -1, {1, 1}, {4}},
                       {2, -1, {2, 0}, {5}},
                       {2, -1, {2, 1}, {6}},
                       {3, -1, {3, 0}, {7}},
                       {3, -1, {3, 1}, {8}}};
  MeshSpecification      out{// The outMesh is distributed and non-contigous
                        {-1, 0, {0, 0}, {0}},
                        {-1, 1, {0, 1}, {0}},
                        {-1, 1, {1, 0}, {0}},
                        {-1, 0, {1, 1}, {0}},
                        {-1, 2, {2, 0}, {0}},
                        {-1, 2, {2, 1}, {0}},
                        {-1, 3, {3, 0}, {0}},
                        {-1, 3, {3, 1}, {0}}};
  ReferenceSpecification ref{// Tests for {0, 1, 0, 0, 0, 0, 0, 0} on the first rank,
                             // {0, 0, 2, 3, 0, 0, 0, 0} on the second, ...
                             {0, {1}},
                             {0, {0}},
                             {0, {0}},
                             {0, {4}},
                             {0, {0}},
                             {0, {0}},
                             {0, {0}},
                             {0, {0}},
                             {1, {0}},
                             {1, {2}},
                             {1, {3}},
                             {1, {0}},
                             {1, {0}},
                             {1, {0}},
                             {1, {0}},
                             {1, {0}},
                             {2, {0}},
                             {2, {0}},
                             {2, {0}},
                             {2, {0}},
                             {2, {5}},
                             {2, {6}},
                             {2, {0}},
                             {2, {0}},
                             {3, {0}},
                             {3, {0}},
                             {3, {0}},
                             {3, {0}},
                             {3, {0}},
                             {3, {0}},
                             {3, {7}},
                             {3, {8}}};

  RadialBasisFctMapping<RadialBasisFctSolver<Multiquadrics>> mapping_on(Mapping::CONSERVATIVE, 2, fct, {{false, false, false}}, Polynomial::ON);
  testDistributed(context, mapping_on, in, out, ref, context.rank * 2);
  RadialBasisFctMapping<RadialBasisFctSolver<Multiquadrics>> mapping_sep(Mapping::CONSERVATIVE, 2, fct, {{false, false, false}}, Polynomial::SEPARATE);
  testDistributed(context, mapping_sep, in, out, ref, context.rank * 2);
  RadialBasisFctMapping<RadialBasisFctSolver<Multiquadrics>> mapping_off(Mapping::CONSERVATIVE, 2, fct, {{false, false, false}}, Polynomial::OFF);
  testDistributed(context, mapping_off, in, out, ref, context.rank * 2);
}

/// Tests a non-contigous owner distributed at the outMesh
BOOST_AUTO_TEST_CASE(testDistributedConservative2DV5Vector)
{
  PRECICE_TEST(""_on(4_ranks).setupIntraComm());
  Multiquadrics fct(5.0);

  MeshSpecification                                          in{// Conservative mapping: The inMesh is local
                       {0, -1, {0, 0}, {1, 4}},
                       {0, -1, {0, 1}, {2, 5}},
                       {1, -1, {1, 0}, {3, 6}},
                       {1, -1, {1, 1}, {4, 7}},
                       {2, -1, {2, 0}, {5, 8}},
                       {2, -1, {2, 1}, {6, 9}},
                       {3, -1, {3, 0}, {7, 10}},
                       {3, -1, {3, 1}, {8, 11}}};
  MeshSpecification                                          out{// The outMesh is distributed and non-contigous
                        {-1, 0, {0, 0}, {0, 0}},
                        {-1, 1, {0, 1}, {0, 0}},
                        {-1, 1, {1, 0}, {0, 0}},
                        {-1, 0, {1, 1}, {0, 0}},
                        {-1, 2, {2, 0}, {0, 0}},
                        {-1, 2, {2, 1}, {0, 0}},
                        {-1, 3, {3, 0}, {0, 0}},
                        {-1, 3, {3, 1}, {0, 0}}};
  ReferenceSpecification                                     ref{// Tests for {0, 1, 0, 0, 0, 0, 0, 0} on the first rank,
                             // {0, 0, 2, 3, 0, 0, 0, 0} on the second, ...
                             {0, {1, 4}},
                             {0, {0, 0}},
                             {0, {0, 0}},
                             {0, {4, 7}},
                             {0, {0, 0}},
                             {0, {0, 0}},
                             {0, {0, 0}},
                             {0, {0, 0}},
                             {1, {0, 0}},
                             {1, {2, 5}},
                             {1, {3, 6}},
                             {1, {0, 0}},
                             {1, {0, 0}},
                             {1, {0, 0}},
                             {1, {0, 0}},
                             {1, {0, 0}},
                             {2, {0, 0}},
                             {2, {0, 0}},
                             {2, {0, 0}},
                             {2, {0, 0}},
                             {2, {5, 8}},
                             {2, {6, 9}},
                             {2, {0, 0}},
                             {2, {0, 0}},
                             {3, {0, 0}},
                             {3, {0, 0}},
                             {3, {0, 0}},
                             {3, {0, 0}},
                             {3, {0, 0}},
                             {3, {0, 0}},
                             {3, {7, 10}},
                             {3, {8, 11}}};
  RadialBasisFctMapping<RadialBasisFctSolver<Multiquadrics>> mapping_on(Mapping::CONSERVATIVE, 2, fct, {{false, false, false}}, Polynomial::ON);
  testDistributed(context, mapping_on, in, out, ref, context.rank * 2);
  RadialBasisFctMapping<RadialBasisFctSolver<Multiquadrics>> mapping_sep(Mapping::CONSERVATIVE, 2, fct, {{false, false, false}}, Polynomial::SEPARATE);
  testDistributed(context, mapping_sep, in, out, ref, context.rank * 2);
  RadialBasisFctMapping<RadialBasisFctSolver<Multiquadrics>> mapping_off(Mapping::CONSERVATIVE, 2, fct, {{false, false, false}}, Polynomial::OFF);
  testDistributed(context, mapping_off, in, out, ref, context.rank * 2);
}

void testTagging(const TestContext &context,
                 MeshSpecification  inMeshSpec,
                 MeshSpecification  outMeshSpec,
                 MeshSpecification  shouldTagFirstRound,
                 MeshSpecification  shouldTagSecondRound,
                 bool               consistent)
{
  int meshDimension  = inMeshSpec.at(0).position.size();
  int valueDimension = inMeshSpec.at(0).value.size();

  mesh::PtrMesh inMesh(new mesh::Mesh("InMesh", meshDimension, testing::nextMeshID()));
  mesh::PtrData inData = inMesh->createData("InData", valueDimension, 0_dataID);
  getDistributedMesh(context, inMeshSpec, inMesh, inData);

  mesh::PtrMesh outMesh(new mesh::Mesh("outMesh", meshDimension, testing::nextMeshID()));
  mesh::PtrData outData = outMesh->createData("OutData", valueDimension, 1_dataID);
  getDistributedMesh(context, outMeshSpec, outMesh, outData);

  Gaussian                                              fct(4.5); // Support radius approx. 1
  Mapping::Constraint                                   constr = consistent ? Mapping::CONSISTENT : Mapping::CONSERVATIVE;
  RadialBasisFctMapping<RadialBasisFctSolver<Gaussian>> mapping(constr, 2, fct, {{false, false, false}}, Polynomial::SEPARATE);
  inMesh->computeBoundingBox();
  outMesh->computeBoundingBox();

  mapping.setMeshes(inMesh, outMesh);
  mapping.tagMeshFirstRound();

  for (const auto &v : inMesh->vertices()) {
    auto pos   = std::find_if(shouldTagFirstRound.begin(), shouldTagFirstRound.end(),
                              [meshDimension, &v](const VertexSpecification &spec) {
                              return std::equal(spec.position.data(), spec.position.data() + meshDimension, v.getCoords().data());
                            });
    bool found = pos != shouldTagFirstRound.end();
    BOOST_TEST(found >= v.isTagged(),
               "FirstRound: Vertex " << v << " is tagged, but should not be.");
    BOOST_TEST(found <= v.isTagged(),
               "FirstRound: Vertex " << v << " is not tagged, but should be.");
  }

  mapping.tagMeshSecondRound();

  for (const auto &v : inMesh->vertices()) {
    auto posFirst    = std::find_if(shouldTagFirstRound.begin(), shouldTagFirstRound.end(),
                                    [meshDimension, &v](const VertexSpecification &spec) {
                                   return std::equal(spec.position.data(), spec.position.data() + meshDimension, v.getCoords().data());
                                 });
    bool foundFirst  = posFirst != shouldTagFirstRound.end();
    auto posSecond   = std::find_if(shouldTagSecondRound.begin(), shouldTagSecondRound.end(),
                                    [meshDimension, &v](const VertexSpecification &spec) {
                                    return std::equal(spec.position.data(), spec.position.data() + meshDimension, v.getCoords().data());
                                  });
    bool foundSecond = posSecond != shouldTagSecondRound.end();
    BOOST_TEST(foundFirst <= v.isTagged(), "SecondRound: Vertex " << v
                                                                  << " is not tagged, but should be from the first round.");
    BOOST_TEST(foundSecond <= v.isTagged(), "SecondRound: Vertex " << v
                                                                   << " is not tagged, but should be.");
    BOOST_TEST((foundSecond or foundFirst) >= v.isTagged(), "SecondRound: Vertex " << v
                                                                                   << " is tagged, but should not be.");
  }
}

BOOST_AUTO_TEST_CASE(testTagFirstRound)
{
  PRECICE_TEST(""_on(4_ranks).setupIntraComm())
  //    *
  //    + <-- owned
  //* * x * *
  //    *
  //    *
  MeshSpecification outMeshSpec = {
      {0, -1, {0, 0}, {0}}};
  MeshSpecification inMeshSpec = {
      {0, -1, {-1, 0}, {1}}, // inside
      {0, -1, {-2, 0}, {1}}, // outside
      {0, 0, {1, 0}, {1}},   // inside, owner
      {0, -1, {2, 0}, {1}},  // outside
      {0, -1, {0, -1}, {1}}, // inside
      {0, -1, {0, -2}, {1}}, // outside
      {0, -1, {0, 1}, {1}},  // inside
      {0, -1, {0, 2}, {1}}   // outside
  };
  MeshSpecification shouldTagFirstRound = {
      {0, -1, {-1, 0}, {1}},
      {0, -1, {1, 0}, {1}},
      {0, -1, {0, -1}, {1}},
      {0, -1, {0, 1}, {1}}};
  MeshSpecification shouldTagSecondRound = {
      {0, -1, {2, 0}, {1}}};
  testTagging(context, inMeshSpec, outMeshSpec, shouldTagFirstRound, shouldTagSecondRound, true);
  // For conservative just swap meshes.
  testTagging(context, outMeshSpec, inMeshSpec, shouldTagFirstRound, shouldTagSecondRound, false);
}

BOOST_AUTO_TEST_SUITE_END() // Parallel

BOOST_AUTO_TEST_SUITE(Serial)

#undef doLocalCode
#define doLocalCode(Type, function, polynomial)                                                                                                                    \
  {                                                                                                                                                                \
    RadialBasisFctMapping<RadialBasisFctSolver<Type>> consistentMap2D(Mapping::CONSISTENT, 2, function, {{false, false, false}}, polynomial);                      \
    perform2DTestConsistentMapping(consistentMap2D);                                                                                                               \
    RadialBasisFctMapping<RadialBasisFctSolver<Type>> consistentMap2DVector(Mapping::CONSISTENT, 2, function, {{false, false, false}}, polynomial);                \
    perform2DTestConsistentMappingVector(consistentMap2DVector);                                                                                                   \
    RadialBasisFctMapping<RadialBasisFctSolver<Type>> consistentMap3D(Mapping::CONSISTENT, 3, function, {{false, false, false}}, polynomial);                      \
    perform3DTestConsistentMapping(consistentMap3D);                                                                                                               \
    RadialBasisFctMapping<RadialBasisFctSolver<Type>> scaledConsistentMap2D(Mapping::SCALED_CONSISTENT_SURFACE, 2, function, {{false, false, false}}, polynomial); \
    perform2DTestScaledConsistentMapping(scaledConsistentMap2D);                                                                                                   \
    RadialBasisFctMapping<RadialBasisFctSolver<Type>> scaledConsistentMap3D(Mapping::SCALED_CONSISTENT_SURFACE, 3, function, {{false, false, false}}, polynomial); \
    perform3DTestScaledConsistentMapping(scaledConsistentMap3D);                                                                                                   \
    RadialBasisFctMapping<RadialBasisFctSolver<Type>> conservativeMap2D(Mapping::CONSERVATIVE, 2, function, {{false, false, false}}, polynomial);                  \
    perform2DTestConservativeMapping(conservativeMap2D);                                                                                                           \
    RadialBasisFctMapping<RadialBasisFctSolver<Type>> conservativeMap2DVector(Mapping::CONSERVATIVE, 2, function, {{false, false, false}}, polynomial);            \
    perform2DTestConservativeMappingVector(conservativeMap2DVector);                                                                                               \
    RadialBasisFctMapping<RadialBasisFctSolver<Type>> conservativeMap3D(Mapping::CONSERVATIVE, 3, function, {{false, false, false}}, polynomial);                  \
    perform3DTestConservativeMapping(conservativeMap3D);                                                                                                           \
  }

BOOST_AUTO_TEST_CASE(MapThinPlateSplines)
{
  PRECICE_TEST(1_rank);
  ThinPlateSplines fct;
  doLocalCode(ThinPlateSplines, fct, Polynomial::ON);
  doLocalCode(ThinPlateSplines, fct, Polynomial::SEPARATE);
}

BOOST_AUTO_TEST_CASE(MapMultiquadrics)
{
  PRECICE_TEST(1_rank);
  Multiquadrics fct(1e-3);
  doLocalCode(Multiquadrics, fct, Polynomial::ON);
  doLocalCode(Multiquadrics, fct, Polynomial::SEPARATE);
}

BOOST_AUTO_TEST_CASE(MapInverseMultiquadrics)
{
  PRECICE_TEST(1_rank);
  InverseMultiquadrics fct(1e-3);
  doLocalCode(InverseMultiquadrics, fct, Polynomial::SEPARATE);
}

BOOST_AUTO_TEST_CASE(MapVolumeSplines)
{
  PRECICE_TEST(1_rank);
  VolumeSplines fct;
  doLocalCode(VolumeSplines, fct, Polynomial::ON);
  doLocalCode(VolumeSplines, fct, Polynomial::SEPARATE);
}

BOOST_AUTO_TEST_CASE(MapGaussian)
{
  PRECICE_TEST(1_rank);
  Gaussian fct(1.0);
  doLocalCode(Gaussian, fct, Polynomial::SEPARATE);
}

BOOST_AUTO_TEST_CASE(MapCompactThinPlateSplinesC2)
{
  PRECICE_TEST(1_rank);
  double                    supportRadius = 1.2;
  CompactThinPlateSplinesC2 fct(supportRadius);
  doLocalCode(CompactThinPlateSplinesC2, fct, Polynomial::SEPARATE);
}

BOOST_AUTO_TEST_CASE(MapCompactPolynomialC0)
{
  PRECICE_TEST(1_rank);
  double              supportRadius = 1.2;
  CompactPolynomialC0 fct(supportRadius);
  doLocalCode(CompactPolynomialC0, fct, Polynomial::SEPARATE);
}

BOOST_AUTO_TEST_CASE(MapCompactPolynomialC2)
{
  PRECICE_TEST(1_rank);
  double              supportRadius = 1.2;
  CompactPolynomialC2 fct(supportRadius);
  doLocalCode(CompactPolynomialC2, fct, Polynomial::SEPARATE);
}

BOOST_AUTO_TEST_CASE(MapCompactPolynomialC4)
{
  PRECICE_TEST(1_rank);
  double              supportRadius = 1.2;
  CompactPolynomialC4 fct(supportRadius);
  doLocalCode(CompactPolynomialC4, fct, Polynomial::SEPARATE);
}

BOOST_AUTO_TEST_CASE(MapCompactPolynomialC6)
{
  PRECICE_TEST(1_rank);
  double              supportRadius = 1.2;
  CompactPolynomialC6 fct(supportRadius);
  doLocalCode(CompactPolynomialC6, fct, Polynomial::SEPARATE);
}
#undef doLocalCode

void testDeadAxis2d(Polynomial polynomial, Mapping::Constraint constraint)
{
  using Eigen::Vector2d;
  int dimensions = 2;

  bool xDead = false;
  bool yDead = true;
  bool zDead = false;

  ThinPlateSplines                                              fct;
  RadialBasisFctMapping<RadialBasisFctSolver<ThinPlateSplines>> mapping(constraint, dimensions, fct,
                                                                        {{xDead, yDead, zDead}}, polynomial);

  // Create mesh to map from
  mesh::PtrMesh inMesh(new mesh::Mesh("InMesh", dimensions, testing::nextMeshID()));
  mesh::PtrData inData   = inMesh->createData("InData", 1, 0_dataID);
  int           inDataID = inData->getID();
  inMesh->createVertex(Vector2d(0.0, 1.0));
  inMesh->createVertex(Vector2d(1.0, 1.0));
  inMesh->createVertex(Vector2d(2.0, 1.0));
  inMesh->createVertex(Vector2d(3.0, 1.0));
  inMesh->allocateDataValues();
  addGlobalIndex(inMesh);
  inMesh->setGlobalNumberOfVertices(inMesh->vertices().size());

  auto &values = inData->values();
  values << 1.0, 2.0, 2.0, 1.0;

  // Create mesh to map to
  mesh::PtrMesh outMesh(new mesh::Mesh("OutMesh", dimensions, testing::nextMeshID()));
  mesh::PtrData outData   = outMesh->createData("OutData", 1, 1_dataID);
  int           outDataID = outData->getID();
  outMesh->createVertex(Vector2d(0, 1.));
  outMesh->createVertex(Vector2d(3, 1.));
  outMesh->createVertex(Vector2d(1.3, 1.));
  outMesh->createVertex(Vector2d(5, 1.));
  outMesh->allocateDataValues();
  addGlobalIndex(outMesh);
  outMesh->setGlobalNumberOfVertices(outMesh->vertices().size());

  // Setup mapping with mapping coordinates and geometry used
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  BOOST_TEST(mapping.hasComputedMapping() == true);
  if (constraint == Mapping::CONSISTENT) {
    if (polynomial == Polynomial::OFF) {
      BOOST_TEST(testing::equals(outData->values()(2), 2.0522549299731567, 1e-7));
    } else if (polynomial == Polynomial::SEPARATE) {
      BOOST_TEST(testing::equals(outData->values()(2), 2.0896514371485777, 1e-7));
    } else {
      BOOST_TEST(testing::equals(outData->values()(2), 2.1180354377884774, 1e-7));
    }
  } else {
    if (polynomial == Polynomial::OFF) {
      BOOST_TEST(testing::equals(outData->values()(1), 1.8471144693068295, 1e-7));
    } else if (polynomial == Polynomial::SEPARATE) {
      BOOST_TEST(testing::equals(outData->values()(1), 1.8236736422730249, 1e-7));
    } else {
      BOOST_TEST(testing::equals(outData->values()(1), 1.7587181970483183, 1e-7));
    }
  }
}

void testDeadAxis3d(Polynomial polynomial, Mapping::Constraint constraint)
{
  using Eigen::Vector3d;
  int dimensions = 3;

  ThinPlateSplines fct;
  bool             xDead = false;
  bool             yDead = true;
  bool             zDead = false;
  using Mapping          = RadialBasisFctMapping<RadialBasisFctSolver<ThinPlateSplines>>;
  Mapping mapping(constraint, dimensions, fct, {{xDead, yDead, zDead}}, polynomial);

  // Create mesh to map from
  mesh::PtrMesh inMesh(new mesh::Mesh("InMesh", dimensions, testing::nextMeshID()));
  mesh::PtrData inData   = inMesh->createData("InData", 1, 0_dataID);
  int           inDataID = inData->getID();
  inMesh->createVertex(Vector3d(0.0, 3.0, 0.0));
  inMesh->createVertex(Vector3d(1.0, 3.0, 0.0));
  inMesh->createVertex(Vector3d(0.0, 3.0, 1.0));
  inMesh->createVertex(Vector3d(1.0, 3.0, 1.0));
  inMesh->allocateDataValues();
  addGlobalIndex(inMesh);
  inMesh->setGlobalNumberOfVertices(inMesh->vertices().size());

  auto &values = inData->values();
  values << 1.0, 2.0, 3.0, 4.0;

  // Create mesh to map to
  mesh::PtrMesh outMesh(new mesh::Mesh("OutMesh", dimensions, testing::nextMeshID()));
  mesh::PtrData outData   = outMesh->createData("OutData", 1, 1_dataID);
  int           outDataID = outData->getID();
  outMesh->createVertex(Vector3d(0.0, 2.9, 0.0));
  outMesh->createVertex(Vector3d(0.8, 2.9, 0.1));
  outMesh->createVertex(Vector3d(0.1, 2.9, 0.9));
  outMesh->createVertex(Vector3d(1.1, 2.9, 1.1));
  outMesh->allocateDataValues();
  addGlobalIndex(outMesh);
  outMesh->setGlobalNumberOfVertices(outMesh->vertices().size());

  // Setup mapping with mapping coordinates and geometry used
  mapping.setMeshes(inMesh, outMesh);
  BOOST_TEST(mapping.hasComputedMapping() == false);

  mapping.computeMapping();
  mapping.map(inDataID, outDataID);
  BOOST_TEST(mapping.hasComputedMapping() == true);

  if (constraint == Mapping::CONSISTENT) {
    if (polynomial == Polynomial::OFF) {
      const double tolerance = 1e-7;
      BOOST_TEST(outData->values()(0) == 1.0);
      BOOST_TEST(testing::equals(outData->values()(1), -0.454450524334, tolerance));
      BOOST_TEST(testing::equals(outData->values()(2), 0.99146426249, tolerance));
      BOOST_TEST(testing::equals(outData->values()(3), 6.98958304876, tolerance));
    } else {
      BOOST_TEST(outData->values()(0) == 1.0);
      BOOST_TEST(outData->values()(1) == 2.0);
      BOOST_TEST(outData->values()(2) == 2.9);
      BOOST_TEST(outData->values()(3) == 4.3);
    }
  } else {
    if (polynomial == Polynomial::OFF) {
      const double tolerance = 1e-6;
      BOOST_TEST(testing::equals(outData->values()(0), 1.17251596926, tolerance));
      BOOST_TEST(testing::equals(outData->values()(1), 4.10368825944, tolerance));
      BOOST_TEST(testing::equals(outData->values()(2), 3.56931954192, tolerance));
      BOOST_TEST(testing::equals(outData->values()(3), 3.40160932341, tolerance));
    } else if (polynomial == Polynomial::ON) {
      const double tolerance = 1e-6;
      BOOST_TEST(testing::equals(outData->values()(0), 0.856701171969, tolerance));
      BOOST_TEST(testing::equals(outData->values()(1), 2.38947124326, tolerance));
      BOOST_TEST(testing::equals(outData->values()(2), 3.34078733786, tolerance));
      BOOST_TEST(testing::equals(outData->values()(3), 3.41304024691, tolerance));
    } else {
      const double tolerance = 1e-6;
      BOOST_TEST(testing::equals(outData->values()(0), 0.380480856704, tolerance));
      BOOST_TEST(testing::equals(outData->values()(1), 2.83529451713, tolerance));
      BOOST_TEST(testing::equals(outData->values()(2), 3.73088270249, tolerance));
      BOOST_TEST(testing::equals(outData->values()(3), 3.05334192368, tolerance));
    }
  }
}

BOOST_AUTO_TEST_CASE(DeadAxis2Consistent)
{
  PRECICE_TEST(1_rank);
  testDeadAxis2d(Polynomial::ON, Mapping::CONSISTENT);
  testDeadAxis2d(Polynomial::OFF, Mapping::CONSISTENT);
  testDeadAxis2d(Polynomial::SEPARATE, Mapping::CONSISTENT);
}

BOOST_AUTO_TEST_CASE(DeadAxis3DConsistent)
{
  PRECICE_TEST(1_rank);
  testDeadAxis3d(Polynomial::ON, Mapping::CONSISTENT);
  testDeadAxis3d(Polynomial::OFF, Mapping::CONSISTENT);
  testDeadAxis3d(Polynomial::SEPARATE, Mapping::CONSISTENT);
}

BOOST_AUTO_TEST_CASE(DeadAxis2Conservative)
{
  PRECICE_TEST(1_rank);
  testDeadAxis2d(Polynomial::ON, Mapping::CONSERVATIVE);
  testDeadAxis2d(Polynomial::OFF, Mapping::CONSERVATIVE);
  testDeadAxis2d(Polynomial::SEPARATE, Mapping::CONSERVATIVE);
}

BOOST_AUTO_TEST_CASE(DeadAxis3DConervative)
{
  PRECICE_TEST(1_rank);
  testDeadAxis3d(Polynomial::ON, Mapping::CONSERVATIVE);
  testDeadAxis3d(Polynomial::OFF, Mapping::CONSERVATIVE);
  testDeadAxis3d(Polynomial::SEPARATE, Mapping::CONSERVATIVE);
}

BOOST_AUTO_TEST_SUITE_END() // Serial

BOOST_AUTO_TEST_SUITE_END() // RadialBasisFunctionMapping
BOOST_AUTO_TEST_SUITE_END()
