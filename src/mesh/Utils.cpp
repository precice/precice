#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <mesh/Edge.hpp>
#include <mesh/Mesh.hpp>
#include <mesh/Utils.hpp>
#include <utils/IntraComm.hpp>
#include "utils/assertion.hpp"

namespace precice::mesh {

/// Given the data and the mesh, this function returns the surface integral. Assumes no overlap exists for the mesh
Eigen::VectorXd integrateSurface(const PtrMesh &mesh, const Eigen::VectorXd &input)
{
  PRECICE_ASSERT(mesh->nVertices() > 0);
  const int       meshDimensions  = mesh->getDimensions();
  const int       valueDimensions = input.size() / mesh->nVertices();
  const auto     &values          = input;
  Eigen::VectorXd integral        = Eigen::VectorXd::Zero(valueDimensions);

  if (meshDimensions == 2) {
    for (const auto &edge : mesh->edges()) {
      int vertex1 = edge.vertex(0).getID() * valueDimensions;
      int vertex2 = edge.vertex(1).getID() * valueDimensions;
      for (int dim = 0; dim < valueDimensions; ++dim) {
        integral(dim) += 0.5 * edge.getLength() * (values(vertex1 + dim) + values(vertex2 + dim));
      }
    }
  } else {
    for (const auto &face : mesh->triangles()) {
      int vertex1 = face.vertex(0).getID() * valueDimensions;
      int vertex2 = face.vertex(1).getID() * valueDimensions;
      int vertex3 = face.vertex(2).getID() * valueDimensions;

      for (int dim = 0; dim < valueDimensions; ++dim) {
        integral(dim) += (face.getArea() / 3.0) * (values(vertex1 + dim) + values(vertex2 + dim) + values(vertex3 + dim));
      }
    }
  }
  return integral;
}

Eigen::VectorXd integrateVolume(const PtrMesh &mesh, const Eigen::VectorXd &input)
{
  PRECICE_ASSERT(mesh->nVertices() > 0);
  const int       meshDimensions  = mesh->getDimensions();
  const int       valueDimensions = input.size() / mesh->nVertices();
  const auto     &values          = input;
  Eigen::VectorXd integral        = Eigen::VectorXd::Zero(valueDimensions);
  if (meshDimensions == 2) {
    for (const auto &face : mesh->triangles()) {
      int vertex1 = face.vertex(0).getID() * valueDimensions;
      int vertex2 = face.vertex(1).getID() * valueDimensions;
      int vertex3 = face.vertex(2).getID() * valueDimensions;
      for (int dim = 0; dim < valueDimensions; ++dim) {
        integral(dim) += (face.getArea() / 3.0) * (values(vertex1 + dim) + values(vertex2 + dim) + values(vertex3 + dim));
      }
    }
  } else {
    for (const auto &tetra : mesh->tetrahedra()) {
      int vertex1 = tetra.vertex(0).getID() * valueDimensions;
      int vertex2 = tetra.vertex(1).getID() * valueDimensions;
      int vertex3 = tetra.vertex(2).getID() * valueDimensions;
      int vertex4 = tetra.vertex(3).getID() * valueDimensions;

      for (int dim = 0; dim < valueDimensions; ++dim) {
        integral(dim) += (tetra.getVolume() / 4.0) * (values(vertex1 + dim) + values(vertex2 + dim) + values(vertex3 + dim) + values(vertex4 + dim));
      }
    }
  }
  return integral;
}

std::size_t countVerticesInBoundingBox(mesh::PtrMesh mesh, const mesh::BoundingBox &bb)
{
  return std::count_if(mesh->vertices().cbegin(), mesh->vertices().cend(), [&bb](const auto &v) { return bb.contains(v); });
}

std::vector<Eigen::VectorXd> findDuplicateVertexCoords(const Mesh &mesh)
{
  // shortcut: exact coordinate match only, vertices closer than a tolerance are not reported
  std::vector<Vertex::RawCoords> coords;
  coords.reserve(mesh.nVertices());
  for (const auto &v : mesh.vertices()) {
    // NaN breaks the strict weak ordering std::sort requires, and such vertices are no meaningful duplicates
    if (std::all_of(v.rawCoords().begin(), v.rawCoords().end(), [](double c) { return std::isfinite(c); })) {
      coords.push_back(v.rawCoords());
    }
  }
  std::sort(coords.begin(), coords.end());

  std::vector<Eigen::VectorXd> duplicates;
  for (auto it = std::adjacent_find(coords.begin(), coords.end()); it != coords.end();
       it      = std::adjacent_find(it, coords.end())) {
    const auto duplicate = *it;
    duplicates.emplace_back(Eigen::Map<const Eigen::VectorXd>(duplicate.data(), mesh.getDimensions()));
    it = std::find_if(it, coords.end(), [&duplicate](const auto &c) { return c != duplicate; });
  }
  return duplicates;
}

} // namespace precice::mesh
