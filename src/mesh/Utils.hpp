#pragma once

#include <algorithm>
#include <array>
#include <mesh/Edge.hpp>
#include <mesh/Mesh.hpp>
#include <utility>

namespace precice {
namespace mesh {

/** return a pointer to the shared vertex of 2 edges
 *
 * If a and b connect the same vertices, then this will simply return one vertex.
 *
 * @param[in] a pointer to Edge a
 * @param[in] b pointer to Edge b
 *
 * @returns a pointer to a Vertex shared by a and b or nullptr otherwise.
 */
inline Vertex *sharedVertex(Edge &a, Edge &b)
{
  Vertex *a0 = &a.vertex(0);
  Vertex *a1 = &a.vertex(1);
  Vertex *b0 = &b.vertex(0);
  Vertex *b1 = &b.vertex(1);
  if (a0 == b0 || a0 == b1) {
    return a0;
  }
  if (a1 == b0 || a1 == b1) {
    return a1;
  }
  return nullptr;
}

/** Calulates the length of an Edge
 *
 * @param[in] e the edge
 *
 * @returns the distance between both vertices of e 
 */
inline double edgeLength(const Edge &e)
{
  return (e.vertex(0).getCoords() - e.vertex(1).getCoords()).norm();
}

/** Calulates the area of an Triangle
 *
 * @param[in] e the triangle
 *
 * @returns surface area of the triangle. If the points are co-linear, the length of the edge.
 */
inline double triangleArea(const Triangle &t)
{
  Eigen::Vector3d vectorA = t.edge(1).vertex(1).getCoords() - t.edge(1).vertex(0).getCoords();
  Eigen::Vector3d vectorB = t.edge(0).vertex(1).getCoords() - t.edge(0).vertex(0).getCoords();
  // Compute cross-product of vector A and vector B
  auto area = 0.5 * vectorA.cross(vectorB).norm();
  // Check if the vectors are co-linear (Pseudo 3D case)
  if (area < 1e-9) {
    // Area is the longest edge.
    std::vector<double> edgeLengths{edgeLength(t.edge(0)), edgeLength(t.edge(1)), edgeLength(t.edge(2))};
    area = *std::max_element(edgeLengths.begin(), edgeLengths.end());
  }
  return area;
}

template <std::size_t n>
struct Chain {
  /// true if the chain is connected or closed and thus valid
  bool connected;
  /// undefined if not connected
  std::array<Vertex *, n> vertices;
  /// undefined if not connected
  std::array<Edge *, n> edges;
};

/** Generates a chain for an array of edges.
 *
 * The resulting verices and edges are undefined if the chain is not connected.
 * If the edges form a chain, then the first edge of the resulting chain is the first edge of the argument.
 * Also, the first and last vertex of the chain will be the vertices of the first edge of the argument.
 *
 * @param[in] edges an array of pointers to edges to chain together
 *
 * @returns A \ref Chain of the input
 */
template <std::size_t n>
Chain<n> asChain(std::array<mesh::Edge *, n> edges)
{
  static_assert(n > 1, "You already know the answer.");
  Chain<n> chain;
  chain.connected = false;

  // the edge 0 is the starting point
  // find connected edges 1 ... n-2
  for (std::size_t i = 1; i < n - 1; ++i) {
    bool found = false;
    for (std::size_t j = i; j < n; ++j) {
      if (edges[i - 1]->connectedTo(*edges[j])) {
        std::swap(edges[i], edges[j]);
        found = true;
        break;
      }
    }
    if (found == false) {
      return chain;
    }
  }
  // the last edge just needs to be checked
  if (!edges[n - 1]->connectedTo(*edges[n - 2]) ||
      !edges[n - 1]->connectedTo(*edges[0])) {
    return chain;
  }
  chain.edges = edges;

  // now find common vertices
  for (std::size_t i = 0; i < n - 1; ++i) {
    chain.vertices[i] = sharedVertex(*edges[i], *edges[i + 1]);
  }
  chain.vertices[n - 1] = sharedVertex(*edges[0], *edges[n - 1]);
  chain.connected       = true;
  return chain;
}

/// Given a mesh and an array of vertexIDS, this function returns an array of pointers to vertices
template <std::size_t n>
std::array<Vertex *, n> vertexPtrsFor(Mesh &mesh, const std::array<int, n> &vertexIDs)
{
  std::array<Vertex *, n> vptrs;
  std::transform(vertexIDs.begin(), vertexIDs.end(), vptrs.begin(),
                 [&mesh](int id) { return &(mesh.vertices()[id]); });
  return vptrs;
}

/// Given a mesh and an array of vertexIDS, this function returns an array of coordinates of the vertices
template <std::size_t n>
std::array<Eigen::VectorXd, n> coordsFor(const Mesh &mesh, const std::array<int, n> &vertexIDs)
{
  std::array<Eigen::VectorXd, n> coords;
  std::transform(vertexIDs.begin(), vertexIDs.end(), coords.begin(),
                 [&mesh](int id) { return mesh.vertices()[id].getCoords(); });
  return coords;
}

/// Given an array of vertex pointers, this function returns an array of coordinates of the vertices
template <std::size_t n>
std::array<Eigen::VectorXd, n> coordsFor(const std::array<Vertex *, n> &vertexPtrs)
{
  std::array<Eigen::VectorXd, n> coords;
  std::transform(vertexPtrs.begin(), vertexPtrs.end(), coords.begin(),
                 [](Vertex *v) { return v->getCoords(); });
  return coords;
}

/// Given the data and the mesh, this function returns the surface integral. Assumes no overlap exists for the mesh
inline Eigen::VectorXd integrate(PtrMesh mesh, PtrData data)
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
        integral(dim) += 0.5 * edgeLength(edge) * (values(vertex1 + dim) + values(vertex2 + dim));
      }
    }
  } else {
    for (const auto &face : mesh->triangles()) {
      int vertex1 = face.vertex(0).getID() * valueDimensions;
      int vertex2 = face.vertex(1).getID() * valueDimensions;
      int vertex3 = face.vertex(2).getID() * valueDimensions;
      for (int dim = 0; dim < valueDimensions; ++dim) {
        integral(dim) += (triangleArea(face) / 3.0) * (values(vertex1 + dim) + values(vertex2 + dim) + values(vertex3 + dim));
      }
    }
  }
  return std::move(integral);
}

/// Given the data and the mesh, this function returns the surface integral
inline Eigen::VectorXd integrateOverlap(PtrMesh mesh, PtrData data)
{
  const int       valueDimensions = data->getDimensions();
  const int       meshDimensions  = mesh->getDimensions();
  const auto &    values          = data->values();
  Eigen::VectorXd integral        = Eigen::VectorXd::Zero(valueDimensions);

  if (meshDimensions == 2) {
    for (const auto &edge : mesh->edges()) {
      int vertex1 = edge.vertex(0).getID() * valueDimensions;
      int vertex2 = edge.vertex(1).getID() * valueDimensions;
      if (edge.vertex(0).isOwner() and edge.vertex(1).isOwner()) {
        for (int dim = 0; dim < valueDimensions; ++dim) {
          integral(dim) += 0.5 * edgeLength(edge) * (values(vertex1 + dim) + values(vertex2 + dim));
        }
      }
    }
  } else {
    for (const auto &face : mesh->triangles()) {
      int vertex1 = face.vertex(0).getID() * valueDimensions;
      int vertex2 = face.vertex(1).getID() * valueDimensions;
      int vertex3 = face.vertex(2).getID() * valueDimensions;
      if (face.vertex(0).isOwner() and face.vertex(1).isOwner() and face.vertex(2).isOwner()) {
        for (int dim = 0; dim < valueDimensions; ++dim) {
          integral(dim) += (triangleArea(face) / 3.0) * (values(vertex1 + dim) + values(vertex2 + dim) + values(vertex3 + dim));
        }
      }
    }
  }
  return std::move(integral);
}

} // namespace mesh
} // namespace precice
