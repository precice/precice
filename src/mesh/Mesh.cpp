#include "Mesh.hpp"
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <algorithm>
#include <array>
#include <boost/container/flat_map.hpp>
#include "Edge.hpp"
#include "Quad.hpp"
#include "RTree.hpp"
#include "Triangle.hpp"
#include "math/math.hpp"
#include "utils/EigenHelperFunctions.hpp"

namespace precice {
namespace mesh {

Mesh::Mesh(
    const std::string &name,
    int                dimensions,
    bool               flipNormals,
    int                id)
    : _name(name),
      _dimensions(dimensions),
      _flipNormals(flipNormals),
      _id(id)
{
  PRECICE_ASSERT((_dimensions == 2) || (_dimensions == 3), _dimensions);
  PRECICE_ASSERT(_name != std::string(""));

  meshChanged.connect([](Mesh &m) { rtree::clear(m); });
  meshDestroyed.connect([](Mesh &m) { rtree::clear(m); });
}

Mesh::~Mesh()
{
  meshDestroyed(*this); // emit signal
}

Mesh::VertexContainer &Mesh::vertices()
{
  return _vertices;
}

const Mesh::VertexContainer &Mesh::vertices() const
{
  return _vertices;
}

Mesh::EdgeContainer &Mesh::edges()
{
  return _edges;
}

const Mesh::EdgeContainer &Mesh::edges() const
{
  return _edges;
}

Mesh::TriangleContainer &Mesh::triangles()
{
  return _triangles;
}

const Mesh::TriangleContainer &Mesh::triangles() const
{
  return _triangles;
}

Mesh::QuadContainer &Mesh::quads()
{
  return _quads;
}

const Mesh::QuadContainer &Mesh::quads() const
{
  return _quads;
}

int Mesh::getDimensions() const
{
  return _dimensions;
}

Edge &Mesh::createEdge(
    Vertex &vertexOne,
    Vertex &vertexTwo)
{
  _edges.emplace_back(vertexOne, vertexTwo, _manageEdgeIDs.getFreeID());
  return _edges.back();
}

Edge &Mesh::createUniqueEdge(
    Vertex &vertexOne,
    Vertex &vertexTwo)
{
  const std::array<int, 2> vids{vertexOne.getID(), vertexTwo.getID()};
  const auto               eend = edges().end();
  auto                     pos  = std::find_if(edges().begin(), eend,
                          [&vids](const Edge &e) -> bool {
                            const std::array<int, 2> eids{e.vertex(0).getID(), e.vertex(1).getID()};
                            return std::is_permutation(vids.begin(), vids.end(), eids.begin());
                          });
  if (pos != eend) {
    return *pos;
  } else {
    return createEdge(vertexOne, vertexTwo);
  }
}

Triangle &Mesh::createTriangle(
    Edge &edgeOne,
    Edge &edgeTwo,
    Edge &edgeThree)
{
  PRECICE_CHECK(
      edgeOne.connectedTo(edgeTwo) &&
          edgeTwo.connectedTo(edgeThree) &&
          edgeThree.connectedTo(edgeOne),
      "Edges are not connected!");
  _triangles.emplace_back(edgeOne, edgeTwo, edgeThree, _manageTriangleIDs.getFreeID());
  return _triangles.back();
}

Quad &Mesh::createQuad(
    Edge &edgeOne,
    Edge &edgeTwo,
    Edge &edgeThree,
    Edge &edgeFour)
{
  _quads.emplace_back(edgeOne, edgeTwo, edgeThree, edgeFour, _manageQuadIDs.getFreeID());
  return _quads.back();
}

PtrData &Mesh::createData(
    const std::string &name,
    int                dimension)
{
  PRECICE_TRACE(name, dimension);
  for (const PtrData data : _data) {
    PRECICE_CHECK(data->getName() != name,
                  "Data \"" << name << "\" cannot be created twice for "
                            << "mesh \"" << _name << "\"!");
  }
  int     id = Data::getDataCount();
  PtrData data(new Data(name, id, dimension));
  _data.push_back(data);
  return _data.back();
}

const Mesh::DataContainer &Mesh::data() const
{
  return _data;
}

const PtrData &Mesh::data(
    int dataID) const
{
  for (const PtrData &data : _data) {
    if (data->getID() == dataID) {
      return data;
    }
  }
  PRECICE_ERROR("Data with ID = " << dataID << " not found in mesh \"" << _name << "\"!");
}

const std::string &Mesh::getName() const
{
  return _name;
}

bool Mesh::isFlipNormals() const
{
  return _flipNormals;
}

void Mesh::setFlipNormals(
    bool flipNormals)
{
  _flipNormals = flipNormals;
}

int Mesh::getID() const
{
  return _id;
}

bool Mesh::isValidVertexID(int vertexID) const
{
  return (0 <= vertexID) && (static_cast<size_t>(vertexID) < vertices().size());
}

bool Mesh::isValidEdgeID(int edgeID) const
{
  return (0 <= edgeID) && (static_cast<size_t>(edgeID) < edges().size());
}

void Mesh::allocateDataValues()
{
  PRECICE_TRACE(_vertices.size());
  for (PtrData data : _data) {
    int total          = _vertices.size() * data->getDimensions();
    int leftToAllocate = total - data->values().size();
    if (leftToAllocate > 0) {
      utils::append(data->values(), (Eigen::VectorXd) Eigen::VectorXd::Zero(leftToAllocate));
    }
    PRECICE_DEBUG("Data " << data->getName() << " now has " << data->values().size() << " values");
  }
}

void Mesh::computeBoundingBox()
{
  PRECICE_TRACE(_name);
  BoundingBox boundingBox(_dimensions,
                          std::make_pair(std::numeric_limits<double>::max(),
                                         std::numeric_limits<double>::lowest()));
  for (const Vertex &vertex : _vertices) {
    for (int d = 0; d < _dimensions; d++) {
      boundingBox[d].first  = std::min(vertex.getCoords()[d], boundingBox[d].first);
      boundingBox[d].second = std::max(vertex.getCoords()[d], boundingBox[d].second);
    }
  }
  for (int d = 0; d < _dimensions; d++) {
    PRECICE_DEBUG("BoundingBox, dim: " << d << ", first: " << boundingBox[d].first << ", second: " << boundingBox[d].second);
  }
  _boundingBox = std::move(boundingBox);
}

void Mesh::computeState()
{
  PRECICE_TRACE(_name);
  PRECICE_ASSERT(_dimensions == 2 || _dimensions == 3, _dimensions);

  // Compute normals only if faces to derive normal information are available
  size_t size2DFaces = _edges.size();
  size_t size3DFaces = _triangles.size() + _quads.size();
  if (_dimensions == 2 && size2DFaces == 0) {
    return;
  }
  if (_dimensions == 3 && size3DFaces == 0) {
    return;
  }

  // Compute (in 2D) edge normals
  if (_dimensions == 2) {
    for (Edge &edge : _edges) {
      Eigen::VectorXd weightednormal = edge.computeNormal(_flipNormals);

      // Accumulate normal in associated vertices
      for (int i = 0; i < 2; i++) {
        Eigen::VectorXd vertexNormal = edge.vertex(i).getNormal();
        vertexNormal += weightednormal;
        edge.vertex(i).setNormal(vertexNormal);
      }
    }
  }

  if (_dimensions == 3) {
    // Compute normals
    for (Triangle &triangle : _triangles) {
      PRECICE_ASSERT(triangle.vertex(0) != triangle.vertex(1),
                     triangle.vertex(0), triangle.getID());
      PRECICE_ASSERT(triangle.vertex(1) != triangle.vertex(2),
                     triangle.vertex(1), triangle.getID());
      PRECICE_ASSERT(triangle.vertex(2) != triangle.vertex(0),
                     triangle.vertex(2), triangle.getID());

      // Compute normals
      Eigen::VectorXd weightednormal = triangle.computeNormal(_flipNormals);

      // Accumulate area-weighted normal in associated vertices and edges
      for (int i = 0; i < 3; i++) {
        triangle.edge(i).setNormal(triangle.edge(i).getNormal() + weightednormal);
        triangle.vertex(i).setNormal(triangle.vertex(i).getNormal() + weightednormal);
      }
    }

    // Compute quad normals
    for (Quad &quad : _quads) {
      PRECICE_ASSERT(quad.vertex(0) != quad.vertex(1), quad.vertex(0).getCoords(), quad.getID());
      PRECICE_ASSERT(quad.vertex(1) != quad.vertex(2), quad.vertex(1).getCoords(), quad.getID());
      PRECICE_ASSERT(quad.vertex(2) != quad.vertex(3), quad.vertex(2).getCoords(), quad.getID());
      PRECICE_ASSERT(quad.vertex(3) != quad.vertex(0), quad.vertex(3).getCoords(), quad.getID());

      // Compute normals (assuming all vertices are on same plane)
      Eigen::VectorXd weightednormal = quad.computeNormal(_flipNormals);
      // Accumulate area-weighted normal in associated vertices and edges
      for (int i = 0; i < 4; i++) {
        quad.edge(i).setNormal(quad.edge(i).getNormal() + weightednormal);
        quad.vertex(i).setNormal(quad.vertex(i).getNormal() + weightednormal);
      }
    }

    // Normalize edge normals (only done in 3D)
    for (Edge &edge : _edges) {
      // there can be cases when an edge has no adjacent triangle though triangles exist in general (e.g. after filtering)
      edge.setNormal(edge.getNormal().normalized());
    }
  }

  for (Vertex &vertex : _vertices) {
    // there can be cases when a vertex has no edge though edges exist in general (e.g. after filtering)
    vertex.setNormal(vertex.getNormal().normalized());
  }
}

void Mesh::clear()
{
  _quads.clear();
  _triangles.clear();
  _edges.clear();
  _vertices.clear();

  _manageTriangleIDs.resetIDs();
  _manageEdgeIDs.resetIDs();
  _manageVertexIDs.resetIDs();

  meshChanged(*this);

  for (mesh::PtrData data : _data) {
    data->values().resize(0);
  }
}

Mesh::VertexDistribution &Mesh::getVertexDistribution()
{
  return _vertexDistribution;
}

const Mesh::VertexDistribution &Mesh::getVertexDistribution() const
{
  return _vertexDistribution;
}

std::vector<int> &Mesh::getVertexOffsets()
{
  return _vertexOffsets;
}

const std::vector<int> &Mesh::getVertexOffsets() const
{
  return _vertexOffsets;
}

void Mesh::setVertexOffsets(std::vector<int> &vertexOffsets)
{
  _vertexOffsets = vertexOffsets;
}

int Mesh::getGlobalNumberOfVertices() const
{
  return _globalNumberOfVertices;
}

void Mesh::setGlobalNumberOfVertices(int num)
{
  _globalNumberOfVertices = num;
}

void Mesh::addMesh(
    Mesh &deltaMesh)
{
  PRECICE_TRACE();
  PRECICE_ASSERT(_dimensions == deltaMesh.getDimensions());

  boost::container::flat_map<int, Vertex *> vertexMap;
  vertexMap.reserve(deltaMesh.vertices().size());
  Eigen::VectorXd coords(_dimensions);
  for (const Vertex &vertex : deltaMesh.vertices()) {
    coords    = vertex.getCoords();
    Vertex &v = createVertex(coords);
    v.setGlobalIndex(vertex.getGlobalIndex());
    if (vertex.isTagged())
      v.tag();
    v.setOwner(vertex.isOwner());
    PRECICE_ASSERT(vertex.getID() >= 0, vertex.getID());
    vertexMap[vertex.getID()] = &v;
  }

  boost::container::flat_map<int, Edge *> edgeMap;
  edgeMap.reserve(deltaMesh.edges().size());
  // you cannot just take the vertices from the edge and add them,
  // since you need the vertices from the new mesh
  // (which may differ in IDs)
  for (const Edge &edge : deltaMesh.edges()) {
    int vertexIndex1 = edge.vertex(0).getID();
    int vertexIndex2 = edge.vertex(1).getID();
    PRECICE_ASSERT((vertexMap.count(vertexIndex1) == 1) &&
                   (vertexMap.count(vertexIndex2) == 1));
    Edge &e               = createEdge(*vertexMap[vertexIndex1], *vertexMap[vertexIndex2]);
    edgeMap[edge.getID()] = &e;
  }

  if (_dimensions == 3) {
    for (const Triangle &triangle : deltaMesh.triangles()) {
      int edgeIndex1 = triangle.edge(0).getID();
      int edgeIndex2 = triangle.edge(1).getID();
      int edgeIndex3 = triangle.edge(2).getID();
      PRECICE_ASSERT((edgeMap.count(edgeIndex1) == 1) &&
                     (edgeMap.count(edgeIndex2) == 1) &&
                     (edgeMap.count(edgeIndex3) == 1));
      createTriangle(*edgeMap[edgeIndex1], *edgeMap[edgeIndex2], *edgeMap[edgeIndex3]);
    }
  }
  meshChanged(*this);
}

const Mesh::BoundingBox Mesh::getBoundingBox() const
{
  return _boundingBox;
}

const std::vector<double> Mesh::getCOG() const
{
  std::vector<double> cog(_dimensions);
  for (int d = 0; d < _dimensions; d++) {
    cog[d] = (_boundingBox[d].second - _boundingBox[d].first) / 2.0 + _boundingBox[d].first;
  }
  return cog;
}

bool Mesh::operator==(const Mesh &other) const
{
  bool equal = true;
  equal &= _vertices.size() == other._vertices.size() &&
           std::is_permutation(_vertices.begin(), _vertices.end(), other._vertices.begin());
  equal &= _edges.size() == other._edges.size() &&
           std::is_permutation(_edges.begin(), _edges.end(), other._edges.begin());
  equal &= _triangles.size() == other._triangles.size() &&
           std::is_permutation(_triangles.begin(), _triangles.end(), other._triangles.begin());
  equal &= _quads.size() == other._quads.size() &&
           std::is_permutation(_quads.begin(), _quads.end(), other._quads.begin());
  return equal;
}

bool Mesh::operator!=(const Mesh &other) const
{
  return !(*this == other);
}

std::ostream &operator<<(std::ostream &os, const Mesh &m)
{
  os << "Mesh \"" << m.getName() << "\", dimensionality = " << m.getDimensions() << ":\n";
  os << "GEOMETRYCOLLECTION(\n";
  const auto  token = ", ";
  const auto *sep   = "";
  for (auto &vertex : m.vertices()) {
    os << sep << vertex;
    sep = token;
  }
  sep = ",\n";
  for (auto &edge : m.edges()) {
    os << sep << edge;
    sep = token;
  }
  sep = ",\n";
  for (auto &triangle : m.triangles()) {
    os << sep << triangle;
    sep = token;
  }
  sep = ",\n";
  for (auto &quad : m.quads()) {
    os << sep << quad;
    sep = token;
  }
  os << "\n)";
  return os;
}

void Mesh::computeQuadConvexityFromPoints(std::array<int,4> &hull, int startID) const
{

  PRECICE_INFO("Entering compute quad: ");
  int vertexOrderIDs[4];    // Initial ordering of vertices from solver adapter

  // All points need to be projected into a new plane with only 2 coordinates, x and y. These are used to check the quad.
  // These are store in vx and vy. For now keep same valiues
  double vx[4];     // x coordinate of projection onto 2D plane
  double vy[4];     // y coordinate of projection onto 2D plane
  double vz[4];     // y coordinate of projection onto 2D plane

  
  // Valid Quad vertices on a 2D plane
  /*
  vx[0] = 0.5;
  vx[1] = 0.62;
  vx[2] = 0.6;
  vx[3] = 0.3;
  vy[0] = 0.34;
  vy[1] = 0.32;
  vy[2] = 0.24;
  vy[3] = 0.22;
  
  // Invalid Quad vertices on a 2D plane
  vx[0] = 0.5;
  vx[1] = 0.62;
  vx[2] = 0.52;
  vx[3] = 0.51;
  vy[0] = 0.34;
  vy[1] = 0.32;
  vy[2] = 0.31;
  vy[3] = 0.22;
  */

  for (int i = 0; i < 4; i++){
    vertexOrderIDs[i] = hull[i];
    vx[i] = vertices()[vertexOrderIDs[i]].getCoords()[0];
    vy[i] = vertices()[vertexOrderIDs[i]].getCoords()[1];
    vz[i] = vertices()[vertexOrderIDs[i]].getCoords()[1];
  }
  PRECICE_INFO("Vertex IDs are: " << hull[0] << " " << hull[1] << " " << hull[2] << " " << hull[3]);
  PRECICE_INFO("X coordinates are: " << vx[0] << " " << vx[1] << " " << vx[2] << " " << vx[3]);
  PRECICE_INFO("Y coordinates are: " << vy[0] << " " << vy[1] << " " << vy[2] << " " << vy[3]);
  
  //First find point with smallest x coord. This point must be in the convex set then and is the starting point of gift wrapping algorithm
  int idLowestPoint = 0;
  for (int i = 1; i < 4; i++) {
    if (vx[i] < vx[idLowestPoint]){
      idLowestPoint = i;
    }
  }

  // Found starting point. Add this as the first vertex in the convex hull.
  // current is the origin point => hull[0]
  int validVertexIDCounter = 0;           // Counts number of times a valid vertex is found
  int currentVertex = idLowestPoint;      // current valid vertex 
  int nextVertex = 0;                     // Next potential valid vertex
  do
  {
    // Add current point to result
    hull[validVertexIDCounter]=currentVertex;
    
    // Search for a point 'nextVertex' such that orientation(p, x,
    // nextVertex) is clockwise for all points 'x'. The idea
    // is to keep track of last visited most clock-
    // wise point in nextVertex. If any point 'i' is more clock-
    // wise than nextVertex, then update nextVertex.
    nextVertex = (currentVertex + 1)%4;              // remainder resets loop through vector of points
    for (int i = 0; i < 4; i++)
    {
      // If i is more clockwise than, nextVertex, then
      // update nextVertex

      double y1 = vy[currentVertex] - vy[nextVertex];
      double y2 = vy[currentVertex] - vy[i];
      double x1 = vx[currentVertex] - vx[nextVertex];
      double x2 = vx[currentVertex] - vx[i];
      double val = y2 * x1 - y1 * x2;
      //PRECICE_INFO("y1: y2: x1 : x2: val = y2*x1 - y1*x2 = " << val << " " << y1 << " " << y2 << " " << x1 << " " << x2);

      if (val > 0){
        nextVertex = i; 	// clock or counterclock wise
          //PRECICE_INFO("Changeing nextVertex");
      }
    }
    // Now nextVertex is the most clockwise with respect to current
    // Set current as nextVertex for next iteration, so that nextVertex is added to
    // result 'hull'
    currentVertex = nextVertex;
    validVertexIDCounter++; 
  } while (currentVertex != idLowestPoint);  // While we don't come to first point

  if (validVertexIDCounter < 4){
    //Error, quad is invalid
    PRECICE_INFO("Invalid Quad. Hull: " << hull);
  } else {
    PRECICE_INFO("Valid Quad. Hull: " << hull);
  }


  //Ordering of quad is hull 0-1-2-3-0

  // Need to check shortest diagonal
  Eigen::VectorXd distance1(_dimensions),distance2(_dimensions);

  distance1 = vertices()[hull[0]].getCoords() - vertices()[hull[2]].getCoords();
  distance2 = vertices()[hull[1]].getCoords() - vertices()[hull[3]].getCoords();
  PRECICE_INFO("Distance 1: " << distance1.norm() << " and Distance 2: " << distance2.norm());
    
  if (distance1.norm() > distance2.norm()){
    startID = 1;
  }
  
}
  

} // namespace mesh
} // namespace precice
