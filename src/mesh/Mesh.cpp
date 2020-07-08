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

bool Mesh::computeQuadConvexityFromPoints(std::array<int,4> &vertexList) const
{
  /*
    All points need to be projected into a new plane with only 2 coordinates, x' and y'. These are used to check 
    the convexity of the quad. These new coordinates are stored in 'coords'. Vertices in 2D does not need to
    be transformed.
  */
  Eigen::Vector3d coords[4];

  // Normal of the plane of first three points in the list of vertices
  Eigen::Vector3d e_1 = vertices()[vertexList[1]].getCoords() - vertices()[vertexList[0]].getCoords();
  Eigen::Vector3d e_2 = vertices()[vertexList[2]].getCoords() - vertices()[vertexList[0]].getCoords();
  Eigen::Vector3d normalVector = e_1.cross(e_2);

  //Transform Coordinates - coord[0] is the origin
  for (int i = 0; i < 4; i++){
    Eigen::Vector3d distance = vertices()[vertexList[i]].getCoords() - vertices()[vertexList[0]].getCoords();
    coords[i][0] = e_1.dot(distance);
    coords[i][1] = e_2.dot(distance);
    coords[i][2] = normalVector.dot(distance);
  }

  PRECICE_DEBUG("Vertex IDs are: " << vertexList[0] << " " << vertexList[1] << " " << vertexList[2] << " " << vertexList[3]);
  PRECICE_DEBUG("X coordinates are: " << coords[0][0] << " " << coords[1][0] << " " << coords[2][0] << " " << coords[3][0]);
  PRECICE_DEBUG("Y coordinates are: " << coords[0][1] << " " << coords[1][1] << " " << coords[2][1] << " " << coords[3][1]);

  /*
  For the convex hull algorithm, the most left hand point regarding the x coordinate is chosen as the starting point. 
  The algorithm moves in an anti-clockwise position, finding the most right hand coordinate from the 
  previous most right hand point. The algorithm must find 3 other points in order to be a valid quad.
  */
  
  //First find point with smallest x coord. This point must be in the convex set then and is the starting point of gift wrapping algorithm
  int idLowestPoint = 0;
  for (int i = 1; i < 4; i++) {
    if (coords[i][0] < coords[idLowestPoint][0]){
      idLowestPoint = i;
    }
  }

  // Found starting point. Add this as the first vertex in the convex hull.
  // current is the origin point => hull[0]
  int validVertexIDCounter = 0;           // Counts number of times a valid vertex is found. Must be 4 for a valid quad.
  int currentVertex = idLowestPoint;      // current valid vertex 
  int nextVertex = 0;                     // Next potential valid vertex
  do
  {
    // Add current point to result
    vertexList[validVertexIDCounter]=currentVertex;
    
    nextVertex = (currentVertex + 1)%4;              // remainder resets loop through vector of points
    for (int i = 0; i < 4; i++)
    {
      double y1 = coords[currentVertex][1] - coords[nextVertex][1];
      double y2 = coords[currentVertex][1] - coords[i][1];
      double x1 = coords[currentVertex][0] - coords[nextVertex][0];
      double x2 = coords[currentVertex][0] - coords[i][0];
      double val = y2 * x1 - y1 * x2;

      if (val > 0){
        nextVertex = i;
      }
    }
    // Now nextVertex is the most anti-clockwise with respect to current
    // Set current as nextVertex for next iteration, so that nextVertex is added to
    // result 'hull'
    currentVertex = nextVertex;
    validVertexIDCounter++; 
  } while (currentVertex != idLowestPoint);  // While we don't come to first point


  if (validVertexIDCounter < 4){
    //Error, quad is invalid
    PRECICE_DEBUG("Invalid Quad. Number of points in convex hull: " << validVertexIDCounter);
  } else {
    PRECICE_DEBUG("Valid Quad. Number of points in convex hull: " << validVertexIDCounter);
  }

  //Ordering of quad is hull 0-1-2-3-0
  return (validVertexIDCounter == 4);
}

std::array<int,4> Mesh::computeQuadEdgeOrder(std::array<int,4> &edgeList) const
{
  /*
  The first edge in the list is treated as the first edge (edge[0]). An edge that does not share a vertex
  with the frist edge is the 3rd edge (edge[2]). The edge that shares the first vertex with the first
  edge (shares vertex[0] with edge[0]) is the 4th edge (edge[3]) as long as it does not share
  the other vertex with edge[0] (then there is a duplicated edge). The edge sharing vertex[1] and not
  vertex[0] is the 2nd edge (edge[1]). This also provides the vertex ordering to determine the diagonal, 
  but not convexity.
  */
  int edgeOrder[4]; // Will contain new order of edges to form a closed quad
  std::array<int,4> vertexList;

  // The first two vertices are the points on the first edge in edgeList. The other edges are built around this
  vertexList[0] = edges()[edgeList[0]].vertex(0).getID();
  vertexList[1] = edges()[edgeList[0]].vertex(1).getID();

  for (int j = 1; j < 4; j++){  // looping thorugh edges 2, 3 and 4 in edgeList

    int ID1 = edges()[edgeList[j]].vertex(0).getID();
    int ID2 = edges()[edgeList[j]].vertex(1).getID();

    if (ID1 != vertexList[0] && ID2 != vertexList[0] && ID1 != vertexList[1] && ID2 != vertexList[1] ){
        // Doesnt match any point on edge 1, therefore must be edge 3 for the reordered set
        edgeOrder[2] = edgeList[j];
      } else if ((ID1 == vertexList[0] || ID2 == vertexList[0]) && (ID1 != vertexList[1] && ID2 != vertexList[1])){
        // One of the vertices matches V0, and does not match V1, must be edge 4 (connected V3 to V0)
        edgeOrder[3] = edgeList[j];
        if (ID1 == vertexList[0]){
          vertexList[3] = edges()[edgeList[j]].vertex(1).getID();
        }else{
          vertexList[3] = edges()[edgeList[j]].vertex(0).getID();
        }
      }else if((ID1 == vertexList[1] || ID2 == vertexList[1]) && (ID1 != vertexList[0] && ID2 != vertexList[0])){
        // Must be edge 2, as it matches V1 and not V0
        edgeOrder[1] = edgeList[j];
        if (ID1 == vertexList[1]){
          vertexList[2] = edges()[edgeList[j]].vertex(1).getID();
        }else{
          vertexList[2] = edges()[edgeList[j]].vertex(0).getID();
        }
      }
    }

    edgeList[1] = edgeOrder[1];
    edgeList[2] = edgeOrder[2];
    edgeList[3] = edgeOrder[3];

    return vertexList;

}
  
} // namespace mesh
} // namespace precice
