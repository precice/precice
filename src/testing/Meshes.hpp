#pragma once

#include "Eigen/Core"
#include "Testing.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"

namespace precice::testing {

inline void addDummyVertices(size_t nVertices, mesh::Mesh &mesh)
{
  for (size_t i = 0; i < nVertices; ++i) {
    if (mesh.getDimensions() == 2) {
      mesh.createVertex(Eigen::Vector2d(i, 0.0));
    } else {
      mesh.createVertex(Eigen::Vector3d(i, 0.0, 0.0));
    }
  }
}

inline auto makeDummy2DMesh(size_t nVertices)
{
  auto mesh = std::make_shared<mesh::Mesh>("DummyMesh2D", 2, testing::nextMeshID());
  addDummyVertices(nVertices, *mesh);
  return mesh;
}

inline auto makeDummy3DMesh(size_t nVertices)
{
  auto mesh = std::make_shared<mesh::Mesh>("DummyMesh3D", 3, testing::nextMeshID());
  addDummyVertices(nVertices, *mesh);
  return mesh;
}

} // namespace precice::testing
