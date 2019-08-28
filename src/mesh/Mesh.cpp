#include "Mesh.hpp"
#include "Edge.hpp"
#include "Triangle.hpp"
#include "Quad.hpp"
#include "PropertyContainer.hpp"
#include "utils/EigenHelperFunctions.hpp"
#include "math/math.hpp"
#include <Eigen/Core>
#include <Eigen/Geometry>
#include "RTree.hpp"
#include <array>
#include <algorithm>

namespace precice {
namespace mesh {

std::unique_ptr<utils::ManageUniqueIDs> Mesh::_managePropertyIDs;

void Mesh:: resetGeometryIDsGlobally()
{
  if (_managePropertyIDs) {
    _managePropertyIDs->resetIDs();
  }
}

Mesh:: Mesh
(
  const std::string& name,
  int                dimensions,
  bool               flipNormals )
:
  _name(name),
  _dimensions(dimensions),
  _flipNormals(flipNormals)
{
  if (not _managePropertyIDs) {
    _managePropertyIDs.reset(new utils::ManageUniqueIDs);
  }
  PRECICE_ASSERT((_dimensions == 2) || (_dimensions == 3), _dimensions);
  PRECICE_ASSERT(_name != std::string(""));
  _nameIDPairs[_name] = _managePropertyIDs->getFreeID ();
  setProperty(INDEX_GEOMETRY_ID, _nameIDPairs[_name]);

  meshChanged.connect([](Mesh & m){rtree::clear(m);});
  meshDestroyed.connect([](Mesh & m){rtree::clear(m);});
}

Mesh:: ~Mesh()
{
  _content.quads().deleteElements();
  _content.triangles().deleteElements();
  _content.edges().deleteElements();
  _content.vertices().deleteElements();
  _propertyContainers.deleteElements();

  meshDestroyed(*this); // emit signal
}

const Group& Mesh:: content() const
{
  return _content;
}

Mesh::VertexContainer& Mesh:: vertices()
{
   return _content.vertices();
}

const Mesh::VertexContainer& Mesh:: vertices() const
{
   return _content.vertices();
}

Mesh::EdgeContainer& Mesh:: edges()
{
   return _content.edges();
}

const Mesh::EdgeContainer& Mesh:: edges() const
{
   return _content.edges();
}

Mesh::TriangleContainer& Mesh:: triangles()
{
   return _content.triangles();
}

const Mesh::TriangleContainer& Mesh:: triangles() const
{
  return _content.triangles();
}

Mesh::PropertyContainerContainer& Mesh:: propertyContainers()
{
  return _propertyContainers;
}

Mesh::QuadContainer& Mesh:: quads()
{
   return _content.quads();
}

const Mesh::QuadContainer& Mesh:: quads() const
{
  return _content.quads();
}

const Mesh::PropertyContainerContainer& Mesh:: propertyContainers() const
{
  return _propertyContainers;
}

int Mesh:: getDimensions() const
{
  return _dimensions;
}

Edge& Mesh:: createEdge
(
  Vertex& vertexOne,
  Vertex& vertexTwo )
{
  Edge* newEdge = new Edge(vertexOne, vertexTwo, _manageEdgeIDs.getFreeID());
  newEdge->addParent(*this);
  _content.add(newEdge);
  return *newEdge;
}

Edge& Mesh::createUniqueEdge
(
    Vertex& vertexOne,
    Vertex& vertexTwo
)
{ 
    const std::array<int, 2> vids{vertexOne.getID(), vertexTwo.getID()};
    const auto eend = edges().end();
    auto pos = std::find_if(edges().begin(), eend,
            [&vids](const Edge& e) -> bool {
                const std::array<int, 2> eids{e.vertex(0).getID(), e.vertex(1).getID()};
                return std::is_permutation(vids.begin(), vids.end(), eids.begin());
            });
    if (pos != eend) {
        return *pos;
    } else {
        return createEdge(vertexOne, vertexTwo);
    }
}

Triangle& Mesh:: createTriangle
(
  Edge& edgeOne,
  Edge& edgeTwo,
  Edge& edgeThree )
{
  PRECICE_CHECK(
          edgeOne.connectedTo(edgeTwo) &&
          edgeTwo.connectedTo(edgeThree) &&
          edgeThree.connectedTo(edgeOne),
          "Edges are not connected!");
  Triangle* newTriangle = new Triangle (
      edgeOne, edgeTwo, edgeThree, _manageTriangleIDs.getFreeID());
  newTriangle->addParent(*this);
  _content.add(newTriangle);
  return *newTriangle;
}

Quad& Mesh:: createQuad
(
  Edge& edgeOne,
  Edge& edgeTwo,
  Edge& edgeThree,
  Edge& edgeFour )
{
  Quad* newQuad = new Quad (
      edgeOne, edgeTwo, edgeThree, edgeFour, _manageQuadIDs.getFreeID());
  newQuad->addParent(*this);
  _content.add(newQuad);
  return *newQuad;
}

PropertyContainer& Mesh:: createPropertyContainer()
{
  PropertyContainer* newPropertyContainer = new PropertyContainer();
  newPropertyContainer->addParent(*this);
  _propertyContainers.push_back(newPropertyContainer);
  return *newPropertyContainer;
}

PtrData& Mesh:: createData
(
  const std::string& name,
  int                dimension )
{
  PRECICE_TRACE(name, dimension);
  for (const PtrData data : _data) {
    PRECICE_CHECK(data->getName() != name,
          "Data \"" << name << "\" cannot be created twice for " << "mesh \"" << _name << "\"!");
  }
  int id = Data::getDataCount();
  PtrData data(new Data(name, id, dimension));
  _data.push_back(data);
  return _data.back();
}

const Mesh::DataContainer& Mesh:: data() const
{
  return _data;
}

const PtrData& Mesh:: data
(
  int dataID ) const
{
  for (const PtrData& data : _data) {
    if (data->getID() == dataID) {
      return data;
    }
  }
  PRECICE_ERROR("Data with ID = " << dataID << " not found in mesh \"" << _name << "\"!" );
}

PropertyContainer& Mesh:: getPropertyContainer
(
  const std::string & subIDName )
{
  PRECICE_TRACE(subIDName);
  PRECICE_ASSERT(_nameIDPairs.count(subIDName) == 1);
  int id = _nameIDPairs[subIDName];
  for (PropertyContainer& cont : _propertyContainers) {
    if (cont.getProperty<int>(cont.INDEX_GEOMETRY_ID) == id){
      return cont;
    }
  }
  PRECICE_ERROR("Unknown sub ID name \"" << subIDName << "\" in mesh \"" << _name << "\"!");
}

const std::string& Mesh:: getName() const
{
  return _name;
}

bool Mesh:: isFlipNormals() const
{
  return _flipNormals;
}

void Mesh:: setFlipNormals
(
  bool flipNormals )
{
  _flipNormals = flipNormals;
}

PropertyContainer& Mesh:: setSubID
(
  const std::string& subIDNamePostfix )
{
  PRECICE_TRACE(subIDNamePostfix);
  PRECICE_CHECK(subIDNamePostfix != std::string(""),
      "Sub ID postfix of mesh \"" << _name << "\" is not allowed to be an empty string!");
  std::string idName(_name + "-" + subIDNamePostfix);
  PRECICE_CHECK(_nameIDPairs.count(idName) == 0,
      "Sub ID postfix of mesh \"" << _name << "\" is already in use!");
  _nameIDPairs[idName] = _managePropertyIDs->getFreeID();
  PropertyContainer * newPropertyContainer = new PropertyContainer();
  newPropertyContainer->setProperty<int>(PropertyContainer::INDEX_GEOMETRY_ID,
    _nameIDPairs[idName]);
  _propertyContainers.push_back(newPropertyContainer);
  return *newPropertyContainer;
}

const std::map<std::string,int>& Mesh:: getNameIDPairs()
{
  return _nameIDPairs;
}

int Mesh:: getID
(
  const std::string& name ) const
{
  PRECICE_ASSERT(_nameIDPairs.count(name) > 0);
  return _nameIDPairs.find(name)->second;
}

int Mesh:: getID() const
{
  std::map<std::string,int>::const_iterator iter = _nameIDPairs.find(_name);
  PRECICE_ASSERT(iter != _nameIDPairs.end());
  return iter->second;
}

bool Mesh::isValidVertexID(int vertexID) const
{
    return (0 <= vertexID) && (static_cast<size_t>(vertexID) < vertices().size());
}

bool Mesh::isValidEdgeID(int edgeID) const
{
    return (0 <= edgeID) && (static_cast<size_t>(edgeID) < edges().size());
}

void Mesh:: allocateDataValues()
{
  PRECICE_TRACE(_content.vertices().size());
  for (PtrData data : _data) {
    int total = _content.vertices().size() * data->getDimensions();
    int leftToAllocate = total - data->values().size();
    if (leftToAllocate > 0){
      utils::append(data->values(), (Eigen::VectorXd) Eigen::VectorXd::Zero(leftToAllocate));
    }
    PRECICE_DEBUG("Data " << data->getName() << " now has " << data->values().size() << " values");
  }
}

void Mesh:: computeNormals()
{
  PRECICE_TRACE(_name);
  // Compute normals only if faces to derive normal information are available
  size_t size2DFaces = _content.edges().size();
  size_t size3DFaces = _content.triangles().size() + _content.quads().size();
  if (_dimensions == 2 && size2DFaces == 0){
      return;
  }
  if (_dimensions == 3 && size3DFaces == 0){
    return;
  }

  // Compute (in 2D) edge normals
  if (_dimensions == 2) {
      for (Edge& edge : _content.edges()) {
          Eigen::VectorXd weightednormal = edge.computeNormal(_flipNormals);

          // Accumulate normal in associated vertices
          for (int i=0; i < 2; i++){
              Eigen::VectorXd vertexNormal = edge.vertex(i).getNormal();
              vertexNormal += weightednormal;
              edge.vertex(i).setNormal(vertexNormal);
          }
      }
  }

  if (_dimensions == 3){
      // Compute normals
      for (Triangle& triangle : _content.triangles()) {
          PRECICE_ASSERT(triangle.vertex(0) != triangle.vertex(1),
                  triangle.vertex(0), triangle.getID());
          PRECICE_ASSERT(triangle.vertex(1) != triangle.vertex(2),
                  triangle.vertex(1), triangle.getID());
          PRECICE_ASSERT(triangle.vertex(2) != triangle.vertex(0),
                  triangle.vertex(2), triangle.getID());

          // Compute normals
          Eigen::VectorXd weightednormal = triangle.computeNormal(_flipNormals);

          // Accumulate area-weighted normal in associated vertices and edges
          for (int i=0; i < 3; i++){
              triangle.edge(i).setNormal(triangle.edge(i).getNormal() + weightednormal);
              triangle.vertex(i).setNormal(triangle.vertex(i).getNormal() + weightednormal);
          }
      }

      // Compute quad normals
      for (Quad& quad : _content.quads()) {
          PRECICE_ASSERT(quad.vertex(0) != quad.vertex(1), quad.vertex(0).getCoords(), quad.getID());
          PRECICE_ASSERT(quad.vertex(1) != quad.vertex(2), quad.vertex(1).getCoords(), quad.getID());
          PRECICE_ASSERT(quad.vertex(2) != quad.vertex(3), quad.vertex(2).getCoords(), quad.getID());
          PRECICE_ASSERT(quad.vertex(3) != quad.vertex(0), quad.vertex(3).getCoords(), quad.getID());

          // Compute normals (assuming all vertices are on same plane)
          Eigen::VectorXd weightednormal = quad.computeNormal(_flipNormals);
          // Accumulate area-weighted normal in associated vertices and edges
          for (int i=0; i < 4; i++){
              quad.edge(i).setNormal(quad.edge(i).getNormal() + weightednormal);
              quad.vertex(i).setNormal(quad.vertex(i).getNormal() + weightednormal);
          }
      }

      // Normalize edge normals (only done in 3D)
      for (Edge& edge : _content.edges()) {
          // there can be cases when an edge has no adjacent triangle though triangles exist in general (e.g. after filtering)
          edge.setNormal(edge.getNormal().normalized());
      }
  }

  for (Vertex& vertex : _content.vertices()) {
      // there can be cases when a vertex has no edge though edges exist in general (e.g. after filtering)
      vertex.setNormal(vertex.getNormal().normalized());
  }
}

void Mesh:: computeBoundingBox()
{
  PRECICE_TRACE(_name);
  BoundingBox boundingBox(_dimensions,
                              std::make_pair(std::numeric_limits<double>::max(),
                                             std::numeric_limits<double>::lowest()));
  for (const Vertex& vertex : _content.vertices()) {
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

void Mesh:: computeState()
{
  PRECICE_TRACE(_name);
  PRECICE_ASSERT(_dimensions==2 || _dimensions==3, _dimensions);

  computeNormals();
  computeBoundingBox();
}

    
void Mesh:: clear()
{
  _content.quads().deleteElements();
  _content.triangles().deleteElements();
  _content.edges().deleteElements();
  _content.vertices().deleteElements();
  _propertyContainers.deleteElements();

  _content.clear();
  _propertyContainers.clear();

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

void Mesh:: addMesh(
    Mesh& deltaMesh)
{
  PRECICE_TRACE();
  PRECICE_ASSERT(_dimensions==deltaMesh.getDimensions());

  std::map<int, Vertex*> vertexMap;
  std::map<int, Edge*> edgeMap;

  Eigen::VectorXd coords(_dimensions);
  for ( const Vertex& vertex : deltaMesh.vertices() ){
    coords = vertex.getCoords();
    Vertex& v = createVertex (coords);
    v.setGlobalIndex(vertex.getGlobalIndex());
    if(vertex.isTagged()) v.tag();
    v.setOwner(vertex.isOwner());
    PRECICE_ASSERT( vertex.getID() >= 0, vertex.getID() );
    vertexMap[vertex.getID()] = &v;
  }

  // you cannot just take the vertices from the edge and add them,
  // since you need the vertices from the new mesh
  // (which may differ in IDs)
  for (const Edge& edge : deltaMesh.edges()) {
    int vertexIndex1 = edge.vertex(0).getID();
    int vertexIndex2 = edge.vertex(1).getID();
    PRECICE_ASSERT( vertexMap.find(vertexIndex1) != vertexMap.end() );
    PRECICE_ASSERT( vertexMap.find(vertexIndex2) != vertexMap.end() );
    Edge& e = createEdge(*vertexMap[vertexIndex1], *vertexMap[vertexIndex2]);
    edgeMap[edge.getID()] = &e;
  }

  if(_dimensions==3){
    for (const Triangle& triangle : deltaMesh.triangles() ) {
      int edgeIndex1 = triangle.edge(0).getID();
      int edgeIndex2 = triangle.edge(1).getID();
      int edgeIndex3 = triangle.edge(2).getID();
      PRECICE_ASSERT( edgeMap.find(edgeIndex1) != edgeMap.end() );
      PRECICE_ASSERT( edgeMap.find(edgeIndex2) != edgeMap.end() );
      PRECICE_ASSERT( edgeMap.find(edgeIndex3) != edgeMap.end() );
      createTriangle(*edgeMap[edgeIndex1],*edgeMap[edgeIndex2],*edgeMap[edgeIndex3]);
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

bool Mesh::operator==(const Mesh& other) const
{
    auto& myContent = _content;
    auto& otherContent = other._content;
    bool equal = true;
    equal &= myContent.vertices().size() == otherContent.vertices().size() &&
        std::is_permutation(myContent.vertices().begin(), myContent.vertices().end(), otherContent.vertices().begin());
    equal &= myContent.edges().size() == otherContent.edges().size() &&
        std::is_permutation(myContent.edges().begin(), myContent.edges().end(), otherContent.edges().begin());
    equal &= myContent.triangles().size() == otherContent.triangles().size() &&
        std::is_permutation(myContent.triangles().begin(), myContent.triangles().end(), otherContent.triangles().begin());
    equal &= myContent.quads().size() == otherContent.quads().size() &&
        std::is_permutation(myContent.quads().begin(), myContent.quads().end(), otherContent.quads().begin());
    return equal;
}

bool Mesh::operator!=(const Mesh& other) const
{
    return !(*this == other);
}

std::ostream& operator<<(std::ostream& os, const Mesh& m)
{
  os << "Mesh \"" << m.getName() << "\", dimensionality = " << m.getDimensions() << ":\n";
  os << "GEOMETRYCOLLECTION(\n";
  const auto token = ", ";
  const auto* sep = "";
  for (auto& vertex : m.content().vertices()){
      os << sep << vertex; 
      sep = token;
  }
  sep = ",\n";
  for (auto& edge : m.content().edges()){
      os << sep << edge;
      sep = token;
  }
  sep = ",\n";
  for (auto& triangle : m.content().triangles()){
      os << sep << triangle;
      sep = token;
  }
  sep = ",\n";
  for (auto& quad : m.content().quads()){
      os << sep << quad;
      sep = token;
  }
  os << "\n)";
  return os;
}

}} // namespace precice, mesh
