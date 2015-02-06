// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "Mesh.hpp"
#include "Edge.hpp"
#include "Triangle.hpp"
#include "Quad.hpp"
#include "PropertyContainer.hpp"
#include "utils/Globals.hpp"

namespace precice {
namespace mesh {

tarch::logging::Log Mesh:: _log("precice::mesh::Mesh");

utils::ManageUniqueIDs* Mesh:: _managerPropertyIDs = NULL;

void Mesh:: resetGeometryIDsGlobally()
{
  if (_managerPropertyIDs != NULL){
    _managerPropertyIDs->resetIDs();
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
  _flipNormals(flipNormals),
  _nameIDPairs(),
  _content(),
  _data(),
  _manageVertexIDs(),
  _manageEdgeIDs(),
  _manageTriangleIDs(),
  _manageQuadIDs(),
  _listeners(),
  _vertexDistribution(),
  _globalNumberOfVertices(-1),
  _boundingBox()
{
  if (_managerPropertyIDs == NULL){
    _managerPropertyIDs = new utils::ManageUniqueIDs;
  }
  assertion1((_dimensions == 2) || (_dimensions == 3), _dimensions);
  assertion(_name != std::string(""));
  _nameIDPairs[_name] = _managerPropertyIDs->getFreeID ();
  setProperty(INDEX_GEOMETRY_ID, _nameIDPairs[_name]);
}

Mesh:: ~Mesh()
{
  _content.quads().deleteElements();
  _content.triangles().deleteElements();
  _content.edges().deleteElements();
  _content.vertices().deleteElements();
}

const Group& Mesh:: content()
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

void Mesh:: addListener
(
  Mesh::MeshListener& listener )
{
  for (const MeshListener* addedListener : _listeners) {
    if (& listener == addedListener){
      return;
    }
  }
  _listeners.push_back(& listener);
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
  preciceTrace2("createData()", name, dimension);
  for (const PtrData data : _data) {
    preciceCheck(data->getName() != name, "createData()",
                 "Data \"" << name << "\" cannot be created twice for "
                 << "mesh \"" << _name << "\"!");
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
   preciceError("data()", "Data with ID = " << dataID << " not found in mesh \""
                << _name << "\"!" );
}

PropertyContainer& Mesh:: getPropertyContainer
(
  const std::string subIDName )
{
  preciceTrace1("getPropertyContainer()", subIDName);
  assertion(_nameIDPairs.count(subIDName) == 1);
  int id = _nameIDPairs[subIDName];
  for (PropertyContainer& cont : _propertyContainers) {
    if (cont.getProperty<int>(cont.INDEX_GEOMETRY_ID) == id){
      return cont;
    }
  }
  preciceError("getPropertyContainer(string)", "Unknown sub ID name \""
               << subIDName << "\" in mesh \"" << _name << "\"!");
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
  preciceTrace1("setSubID()", subIDNamePostfix);
  preciceCheck(subIDNamePostfix != std::string(""), "setSubID",
      "Sub ID postfix of mesh \"" << _name
      << "\" is not allowed to be an empty string!");
  std::string idName(_name + "-" + subIDNamePostfix);
  preciceCheck(_nameIDPairs.count(idName) == 0, "setSubID",
      "Sub ID postfix of mesh \"" << _name << "\" is already in use!");
  _nameIDPairs[idName] = _managerPropertyIDs->getFreeID();
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
  preciceTrace1("allocateDataValues()", _content.vertices().size());
  for (PtrData data : _data) {
    int total = _content.vertices().size() * data->getDimensions();
    int leftToAllocate = total - data->values().size();
    if (leftToAllocate > 0){
      data->values().append(utils::DynVector(leftToAllocate, 0.0));
    }
    preciceDebug("Data " << data->getName() << " no has "
                 << data->values().size() << " values");
  }
}

void Mesh:: computeState()
{
  preciceTrace("computeState()");
  assertion1(_dimensions==2 || _dimensions==3, _dimensions);
  using utils::DynVector;
  using utils::Vector2D;
  using utils::Vector3D;

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
    assertion1(_dimensions == 3, _dimensions);
    computeNormals = false;
  }

  // Compute edge centers, enclosing radius, and (in 2D) edge normals
  DynVector center(_dimensions);
  DynVector distanceToCenter(_dimensions);
  for (Edge& edge : _content.edges()) {
    center = edge.vertex(0).getCoords();
    center += edge.vertex(1).getCoords();
    center *= 0.5;
    edge.setCenter(center);
    distanceToCenter = edge.vertex(0).getCoords();
    distanceToCenter -= edge.getCenter();
    edge.setEnclosingRadius(tarch::la::norm2(distanceToCenter));
    if (_dimensions == 2 && computeNormals){
      // Compute normal
      Vector2D vectorA = edge.vertex(1).getCoords();
      vectorA -= edge.vertex(0).getCoords();
      Vector2D normal(-1.0 *vectorA[1], vectorA[0]);
      if (not _flipNormals){
        normal *= -1.0; // Invert direction if counterclockwise
      }
      double length = tarch::la::norm2(normal);
      assertion(tarch::la::greater(length, 0.0));
      normal /= length;   // Scale normal vector to length 1
      edge.setNormal(normal);

      // Accumulate normal in associated vertices
      normal *= edge.getEnclosingRadius() * 2.0; // Weight by length
      for (int i=0; i < 2; i++){
        Vector2D vertexNormal = edge.vertex(i).getNormal();
        vertexNormal += normal;
        edge.vertex(i).setNormal(vertexNormal);
      }
    }
  }


  if (_dimensions == 3){
    // Compute triangle centers, radius, and normals
    for (Triangle& triangle : _content.triangles()) {
      assertion2(not tarch::la::equals(triangle.vertex(0).getCoords(),
                 triangle.vertex(1).getCoords()), triangle.vertex(0).getCoords(),
                 triangle.getID());
      assertion2(not tarch::la::equals(triangle.vertex(1).getCoords(),
                 triangle.vertex(2).getCoords()), triangle.vertex(1).getCoords(),
                 triangle.getID());
      assertion2(not tarch::la::equals(triangle.vertex(2).getCoords(),
                 triangle.vertex(0).getCoords()), triangle.vertex(2).getCoords(),
                 triangle.getID());

      // Compute barycenter by using edge centers, since vertex order is not
      // guaranteed.
      Vector3D center;
      center = triangle.edge(0).getCenter();
      center += triangle.edge(1).getCenter();
      center += triangle.edge(2).getCenter();
      center /= 3.0;
      triangle.setCenter(center);

      // Compute enclosing radius centered at barycenter
      Vector3D toCenter;
      toCenter = triangle.getCenter();
      toCenter -= triangle.vertex(0).getCoords();
      double distance0 = tarch::la::norm2(toCenter);
      toCenter = triangle.getCenter();
      toCenter -= triangle.vertex(1).getCoords();
      double distance1 = tarch::la::norm2(toCenter);
      toCenter = triangle.getCenter();
      toCenter -= triangle.vertex(2).getCoords();
      double distance2 = tarch::la::norm2(toCenter);
      double maxDistance = distance0;
      maxDistance = distance1 > maxDistance ? distance1 : maxDistance;
      maxDistance = distance2 > maxDistance ? distance2 : maxDistance;
      triangle.setEnclosingRadius(maxDistance);

      // Compute normals
      if (computeNormals){
        Vector3D vectorA;
        Vector3D vectorB;
        vectorA = triangle.edge(1).getCenter(); // edge() is faster than vertex()
        vectorA -= triangle.edge(0).getCenter();
        vectorB = triangle.edge(2).getCenter();
        vectorB -= triangle.edge(0).getCenter();
        // Compute cross-product of vector A and vector B
        Vector3D normal;
        normal = tarch::la::cross(vectorA, vectorB, normal);
        if ( _flipNormals ){
          normal *= -1.0; // Invert direction if counterclockwise
        }

        // Accumulate area-weighted normal in associated vertices and edges
        for (int i=0; i < 3; i++){
          triangle.edge(i).setNormal(triangle.edge(i).getNormal() + normal);
          triangle.vertex(i).setNormal(triangle.vertex(i).getNormal() + normal);
        }

        // Normalize triangle normal
        double length = tarch::la::norm2(normal);
        normal /= length;
        triangle.setNormal(normal);
      }
    }

    // Compute quad centers, radius, and normals
    for (Quad& quad : _content.quads()) {
      assertion2(not tarch::la::equals(quad.vertex(0).getCoords(),
                 quad.vertex(1).getCoords()), quad.vertex(0).getCoords(),
                 quad.getID());
      assertion2(not tarch::la::equals(quad.vertex(1).getCoords(),
                 quad.vertex(2).getCoords()), quad.vertex(1).getCoords(),
                 quad.getID());
      assertion2(not tarch::la::equals(quad.vertex(2).getCoords(),
                 quad.vertex(3).getCoords()), quad.vertex(2).getCoords(),
                 quad.getID());
      assertion2(not tarch::la::equals(quad.vertex(3).getCoords(),
                 quad.vertex(0).getCoords()), quad.vertex(3).getCoords(),
                 quad.getID());

      // Compute barycenter by using edge centers, since vertex order is not
      // guaranteed.
      Vector3D center;
      center = quad.edge(0).getCenter();
      center += quad.edge(1).getCenter();
      center += quad.edge(2).getCenter();
      center += quad.edge(3).getCenter();
      center /= 4.0;
      quad.setCenter(center);

      // Compute enclosing radius centered at barycenter
      Vector3D toCenter;
      toCenter = quad.getCenter();
      toCenter -= quad.vertex(0).getCoords();
      double distance0 = tarch::la::norm2(toCenter);
      toCenter = quad.getCenter();
      toCenter -= quad.vertex(1).getCoords();
      double distance1 = tarch::la::norm2(toCenter);
      toCenter = quad.getCenter();
      toCenter -= quad.vertex(2).getCoords();
      double distance2 = tarch::la::norm2(toCenter);
      toCenter = quad.getCenter();
      toCenter -= quad.vertex(3).getCoords();
      double distance3 = tarch::la::norm2(toCenter);
      double maxDistance = distance0;
      maxDistance = distance1 > maxDistance ? distance1 : maxDistance;
      maxDistance = distance2 > maxDistance ? distance2 : maxDistance;
      maxDistance = distance3 > maxDistance ? distance3 : maxDistance;
      quad.setEnclosingRadius(maxDistance);


      // Compute normals (assuming all vertices are on same plane)
      if (computeNormals){
        // Two triangles are thought by splitting the quad from vertex 0 to 2.
        // The cross prodcut of the outer edges of the triangles is used to compute
        // the normal direction and area of the triangles. The direction must be
        // the same, while the areas differ in general. The normals are added up
        // and divided by 2 to get the area of the overall quad, since the length
        // does correspond to the parallelogram spanned by the vectors of the
        // cross product, which is twice the area of the corresponding triangles.
        Vector3D vectorA;
        Vector3D vectorB;
        vectorA = quad.vertex(2).getCoords();
        vectorA -= quad.vertex(1).getCoords();
        vectorB = quad.vertex(0).getCoords();
        vectorB -= quad.vertex(1).getCoords();
        // Compute cross-product of vector A and vector B
        Vector3D normal;
        tarch::la::cross(vectorA, vectorB, normal);

        vectorA = quad.vertex(0).getCoords();
        vectorA -= quad.vertex(3).getCoords();
        vectorB = quad.vertex(2).getCoords();
        vectorB -= quad.vertex(3).getCoords();
        Vector3D normalSecondPart;
        tarch::la::cross(vectorA, vectorB, normalSecondPart);

        assertion2(tarch::la::equals(
                   normal/tarch::la::norm2(normal),
                   normalSecondPart/tarch::la::norm2(normalSecondPart)),
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
        double length = tarch::la::norm2(normal);
        normal /= length;
        quad.setNormal(normal);
      }
    }

    // Normalize edge normals (only done in 3D)
    if (computeNormals){
      for (Edge& edge : _content.edges()) {
        double length = tarch::la::norm2(edge.getNormal());
        assertionMsg(tarch::la::greater(length,0.0),
          "Edge vertex coords: (" << edge.vertex(0).getCoords() << "), ("
          << edge.vertex(1).getCoords()
          << "). Hint: Could be inconsistent triangle/quad orientation or "
          << "dangling edge. ");
        edge.setNormal(edge.getNormal() / length);
      }
    }
  }

  // Normalize vertex normals & compute bounding box
  _boundingBox = BoundingBox (_dimensions,
                              std::make_pair(std::numeric_limits<double>::max(),
                                             std::numeric_limits<double>::lowest()));

  for (Vertex& vertex : _content.vertices()) {
    if (computeNormals) {
      double length = tarch::la::norm2(vertex.getNormal());
      // i (benjamin) changed this since there can be cases where a node has no edge though
      // the mesh has edges in general, e.g. after filtering
      //assertion(tarch::la::greater(length,0.0));
      if(tarch::la::greater(length,0.0)){
        vertex.setNormal(vertex.getNormal() / length);
      }
    }
    
    for (int d = 0; d < _dimensions; d++) {
      _boundingBox[d].first  = std::min(vertex.getCoords()[d], _boundingBox[d].first);
      _boundingBox[d].second = std::max(vertex.getCoords()[d], _boundingBox[d].second);
    }
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

  for (mesh::PtrData data : _data) {
    data->values().clear();
  }
}

void Mesh:: notifyListeners()
{
  preciceTrace("notifyListeners()");
  for (MeshListener* listener : _listeners) {
    assertion(listener != NULL);
    listener->meshChanged(*this);
  }
}

void Mesh:: addMesh(
    Mesh& deltaMesh)
{
  preciceTrace("addMesh()");
  assertion(_dimensions==deltaMesh.getDimensions());

  std::map<int, Vertex*> vertexMap;
  std::map<int, Edge*> edgeMap;

  utils::DynVector coords(_dimensions);
  for ( const Vertex& vertex : deltaMesh.vertices() ){
    coords = vertex.getCoords();
    Vertex& v = createVertex (coords);
    assertion1 ( vertex.getID() >= 0, vertex.getID() );
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

}} // namespace precice, mesh
