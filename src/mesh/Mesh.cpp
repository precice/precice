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
  assertion((_dimensions == 2) || (_dimensions == 3), _dimensions);
  assertion(_name != std::string(""));
  _nameIDPairs[_name] = _managePropertyIDs->getFreeID ();
  setProperty(INDEX_GEOMETRY_ID, _nameIDPairs[_name]);

  meshChanged.connect(&rtree::clear);
  meshDestroyed.connect(&rtree::clear);
}

Mesh:: ~Mesh()
{
  _content.quads().deleteElements();
  _content.triangles().deleteElements();
  _content.edges().deleteElements();
  _content.vertices().deleteElements();

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

Triangle& Mesh:: createTriangle
(
  Edge& edgeOne,
  Edge& edgeTwo,
  Edge& edgeThree )
{
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
  TRACE(name, dimension);
  for (const PtrData data : _data) {
    CHECK(data->getName() != name,
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
  ERROR("Data with ID = " << dataID << " not found in mesh \"" << _name << "\"!" );
}

PropertyContainer& Mesh:: getPropertyContainer
(
  const std::string & subIDName )
{
  TRACE(subIDName);
  assertion(_nameIDPairs.count(subIDName) == 1);
  int id = _nameIDPairs[subIDName];
  for (PropertyContainer& cont : _propertyContainers) {
    if (cont.getProperty<int>(cont.INDEX_GEOMETRY_ID) == id){
      return cont;
    }
  }
  ERROR("Unknown sub ID name \"" << subIDName << "\" in mesh \"" << _name << "\"!");
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
  TRACE(subIDNamePostfix);
  CHECK(subIDNamePostfix != std::string(""),
      "Sub ID postfix of mesh \"" << _name << "\" is not allowed to be an empty string!");
  std::string idName(_name + "-" + subIDNamePostfix);
  CHECK(_nameIDPairs.count(idName) == 0,
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
  assertion(_nameIDPairs.count(name) > 0);
  return _nameIDPairs.find(name)->second;
}

int Mesh:: getID() const
{
  std::map<std::string,int>::const_iterator iter = _nameIDPairs.find(_name);
  assertion(iter != _nameIDPairs.end());
  return iter->second;
}

void Mesh:: allocateDataValues()
{
  TRACE(_content.vertices().size());
  for (PtrData data : _data) {
    int total = _content.vertices().size() * data->getDimensions();
    int leftToAllocate = total - data->values().size();
    if (leftToAllocate > 0){
      utils::append(data->values(), (Eigen::VectorXd) Eigen::VectorXd::Zero(leftToAllocate));
    }
    DEBUG("Data " << data->getName() << " no has " << data->values().size() << " values");
  }
}

void Mesh:: computeState()
{
  TRACE(_name);
  assertion(_dimensions==2 || _dimensions==3, _dimensions);

  // Compute normals only if faces to derive normal information are available
  bool computeNormals = true;
  size_t size2DFaces = _content.edges().size();
  size_t size3DFaces = _content.triangles().size() + _content.quads().size();
  if (_dimensions == 2){
    if (size2DFaces == 0){
      computeNormals = false;
    }
  }
  else if (size3DFaces == 0){
    assertion(_dimensions == 3, _dimensions);
    computeNormals = false;
  }

  // Compute edge centers, enclosing radius, and (in 2D) edge normals
  for (Edge& edge : _content.edges()) {
    if (_dimensions == 2 && computeNormals){
      // Compute normal
      Eigen::VectorXd vectorA = edge.vertex(1).getCoords();
      vectorA -= edge.vertex(0).getCoords();
      Eigen::Vector2d normal(-1.0 *vectorA[1], vectorA[0]);
      if (not _flipNormals){
        normal *= -1.0; // Invert direction if counterclockwise
      }
      double length = normal.norm();
      assertion(math::greater(length, 0.0));
      normal /= length;   // Scale normal vector to length 1
      edge.setNormal(normal);

      // Accumulate normal in associated vertices
      normal *= edge.getEnclosingRadius() * 2.0; // Weight by length
      for (int i=0; i < 2; i++){
        Eigen::VectorXd vertexNormal = edge.vertex(i).getNormal();
        vertexNormal += normal;
        edge.vertex(i).setNormal(vertexNormal);
      }
    }
  }

  if (_dimensions == 3){
    // Compute triangle centers, radius, and normals
    for (Triangle& triangle : _content.triangles()) {
      assertion(not math::equals(triangle.vertex(0).getCoords(), triangle.vertex(1).getCoords()),
                triangle.vertex(0).getCoords(),
                triangle.getID());
      assertion(not math::equals(triangle.vertex(1).getCoords(), triangle.vertex(2).getCoords()), triangle.vertex(1).getCoords(),
                triangle.getID());
      assertion(not math::equals(triangle.vertex(2).getCoords(), triangle.vertex(0).getCoords()),
                triangle.vertex(2).getCoords(),
                triangle.getID());

      // Compute normals
      if (computeNormals){
        Eigen::Vector3d vectorA = triangle.edge(1).getCenter() - triangle.edge(0).getCenter(); // edge() is faster than vertex()
        Eigen::Vector3d vectorB = triangle.edge(2).getCenter() - triangle.edge(0).getCenter();
        // Compute cross-product of vector A and vector B
        auto normal = vectorA.cross(vectorB);
        if ( _flipNormals ){
          normal *= -1.0; // Invert direction if counterclockwise
        }

        // Accumulate area-weighted normal in associated vertices and edges
        for (int i=0; i < 3; i++){
          triangle.edge(i).setNormal(triangle.edge(i).getNormal() + normal);
          triangle.vertex(i).setNormal(triangle.vertex(i).getNormal() + normal);
        }

        // Normalize triangle normal
        double length = normal.norm();
        normal /= length;
        triangle.setNormal(normal);
      }
    }

    // Compute quad centers, radius, and normals
    for (Quad& quad : _content.quads()) {
      assertion(not math::equals(quad.vertex(0).getCoords(), quad.vertex(1).getCoords()),
                quad.vertex(0).getCoords(),
                quad.getID());
      assertion(not math::equals(quad.vertex(1).getCoords(), quad.vertex(2).getCoords()),
                quad.vertex(1).getCoords(),
                quad.getID());
      assertion(not math::equals(quad.vertex(2).getCoords(), quad.vertex(3).getCoords()),
                quad.vertex(2).getCoords(),
                quad.getID());
      assertion(not math::equals(quad.vertex(3).getCoords(), quad.vertex(0).getCoords()),
                quad.vertex(3).getCoords(),
                quad.getID());

      // Compute normals (assuming all vertices are on same plane)
      if (computeNormals){
        // Two triangles are thought by splitting the quad from vertex 0 to 2.
        // The cross prodcut of the outer edges of the triangles is used to compute
        // the normal direction and area of the triangles. The direction must be
        // the same, while the areas differ in general. The normals are added up
        // and divided by 2 to get the area of the overall quad, since the length
        // does correspond to the parallelogram spanned by the vectors of the
        // cross product, which is twice the area of the corresponding triangles.
        Eigen::Vector3d vectorA = quad.vertex(2).getCoords() - quad.vertex(1).getCoords();
        Eigen::Vector3d vectorB = quad.vertex(0).getCoords() - quad.vertex(1).getCoords();
        // Compute cross-product of vector A and vector B
        auto normal = vectorA.cross(vectorB);
        
        vectorA = quad.vertex(0).getCoords() - quad.vertex(3).getCoords();
        vectorB = quad.vertex(2).getCoords() - quad.vertex(3).getCoords();
        auto normalSecondPart = vectorA.cross(vectorB);
        
        assertion(math::equals(normal.normalized(), normalSecondPart.normalized()),
                  normal, normalSecondPart);
        normal += normalSecondPart;
        normal *= 0.5;

        if ( _flipNormals ){
          normal *= -1.0; // Invert direction if counterclockwise
        }

        // Accumulate area-weighted normal in associated vertices and edges
        for (int i=0; i < 4; i++){
          quad.edge(i).setNormal(quad.edge(i).getNormal() + normal);
          quad.vertex(i).setNormal(quad.vertex(i).getNormal() + normal);
        }

        // Normalize quad normal
        normal = normal.normalized();
        quad.setNormal(normal);
      }
    }

    // Normalize edge normals (only done in 3D)
    if (computeNormals){
      for (Edge& edge : _content.edges()) {
        double length = edge.getNormal().norm();
        // there can be cases when an edge has no adjacent triangle though triangles exist in general (e.g. after filtering)
        if(math::greater(length,0.0)){
          edge.setNormal(edge.getNormal() / length);
        }
      }
    }
  }

  // Normalize vertex normals & compute bounding box
  _boundingBox = BoundingBox (_dimensions,
                              std::make_pair(std::numeric_limits<double>::max(),
                                             std::numeric_limits<double>::lowest()));

  for (Vertex& vertex : _content.vertices()) {
    if (computeNormals) {
      double length = vertex.getNormal().norm();
      // there can be cases when a vertex has no edge though edges exist in general (e.g. after filtering)
      if(math::greater(length,0.0)){
        vertex.setNormal(vertex.getNormal() / length);
      }
    }
    
    for (int d = 0; d < _dimensions; d++) {
      _boundingBox[d].first  = std::min(vertex.getCoords()[d], _boundingBox[d].first);
      _boundingBox[d].second = std::max(vertex.getCoords()[d], _boundingBox[d].second);
    }
  }
  for (int d = 0; d < _dimensions; d++) {
    DEBUG("BoundingBox, dim: " << d << ", first: " << _boundingBox[d].first << ", second: " << _boundingBox[d].second);
  }
}

    
void Mesh:: clear()
{
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


void Mesh:: addMesh(
    Mesh& deltaMesh)
{
  TRACE();
  assertion(_dimensions==deltaMesh.getDimensions());

  std::map<int, Vertex*> vertexMap;
  std::map<int, Edge*> edgeMap;

  Eigen::VectorXd coords(_dimensions);
  for ( const Vertex& vertex : deltaMesh.vertices() ){
    coords = vertex.getCoords();
    Vertex& v = createVertex (coords);
    v.setGlobalIndex(vertex.getGlobalIndex());
    if(vertex.isTagged()) v.tag();
    v.setOwner(vertex.isOwner());
    assertion ( vertex.getID() >= 0, vertex.getID() );
    vertexMap[vertex.getID()] = &v;
  }

  // you cannot just take the vertices from the edge and add them,
  // since you need the vertices from the new mesh
  // (which may differ in IDs)
  for (const Edge& edge : deltaMesh.edges()) {
    int vertexIndex1 = edge.vertex(0).getID();
    int vertexIndex2 = edge.vertex(1).getID();
    assertion ( vertexMap.find(vertexIndex1) != vertexMap.end() );
    assertion ( vertexMap.find(vertexIndex2) != vertexMap.end() );
    Edge& e = createEdge(*vertexMap[vertexIndex1], *vertexMap[vertexIndex2]);
    edgeMap[edge.getID()] = &e;
  }

  if(_dimensions==3){
    for (const Triangle& triangle : deltaMesh.triangles() ) {
      int edgeIndex1 = triangle.edge(0).getID();
      int edgeIndex2 = triangle.edge(1).getID();
      int edgeIndex3 = triangle.edge(2).getID();
      assertion ( edgeMap.find(edgeIndex1) != edgeMap.end() );
      assertion ( edgeMap.find(edgeIndex2) != edgeMap.end() );
      assertion ( edgeMap.find(edgeIndex3) != edgeMap.end() );
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
  os << "Mesh " << m.getName() << " consisting the Vertices:\n";
  for (auto& vertex : m.content().edges())
      os << "\t" << vertex;
  os << "And the Edges:\n";
  for (auto& edge : m.content().edges())
      os << "\t" << edge;
  os << "And the Triangles:\n";
  for (auto& triangle : m.content().edges())
      os << "\t" << triangle;
  os << "And the Quads:\n";
  for (auto& quad : m.content().quads())
      os << "\t" << quad;
  return os;
}

}} // namespace precice, mesh
