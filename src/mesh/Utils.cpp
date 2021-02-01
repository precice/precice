#include <Eigen/Core>
#include <mesh/Edge.hpp>
#include <mesh/Mesh.hpp>
#include <utils/MasterSlave.hpp>

namespace precice {
namespace mesh {

/// Given the data and the mesh, this function returns the surface integral. Assumes no overlap exists for the mesh
Eigen::VectorXd integrate(PtrMesh mesh, PtrData data)
{
  const int       valueDimensions = data->getDimensions();
  const int       meshDimensions  = mesh->getDimensions();
  const auto &    values          = data->values();
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

} // namespace mesh
} // namespace precice
