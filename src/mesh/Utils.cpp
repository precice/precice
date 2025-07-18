#include <Eigen/Core>
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

std::size_t countVerticesInBoundingBox(mesh::PtrMesh mesh, const std::vector<mesh::BoundingBox> &bbs)
{
  return std::count_if(mesh->vertices().cbegin(),
                       mesh->vertices().cend(),
                       [&bbs](const auto &v) {
                         return std::any_of(bbs.cbegin(),
                                            bbs.cend(),
                                            [&v](const auto &bb) { return bb.contains(v); });
                       });
}

} // namespace precice::mesh
