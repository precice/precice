// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_MESH_MESH_HPP_
#define PRECICE_MESH_MESH_HPP_

#include "mesh/Group.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/PropertyContainer.hpp"
#include "mesh/Data.hpp"
#include "mesh/Vertex.hpp"
#include "tarch/logging/Log.h"
#include "utils/Dimensions.hpp"
#include "utils/PointerVector.hpp"
#include "utils/ManageUniqueIDs.hpp"
#include "boost/utility.hpp"
#include "boost/noncopyable.hpp"
#include "tarch/la/DynamicVector.h"
#include <map>
#include <vector>

namespace precice {
  namespace mesh {
    class PropertyContainer;
  }
}

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace mesh {


/**
 * @brief Container and creator for meshes.
 *
 * A Mesh can consist of Vertices, Edges, Triangles, nested Meshes, and
 * PropertyContainers. It provides functionality to conveniently create those
 * objects.
 *
 * In addition to creating the topological information of a mesh, the Mesh class
 * can also be used to add data to the vertices of a mesh.
 *
 * Usage example: precice::mesh::tests::MeshTest::testDemonstration()
 */
class Mesh : public PropertyContainer, private boost::noncopyable
{
public:

  /**
   * @brief Interface for classes depending on the mesh.
   */
  class MeshListener {
  public:
    virtual void meshChanged ( Mesh& mesh ) =0;
  };

  typedef utils::ptr_vector<Vertex>              VertexContainer;
  typedef utils::ptr_vector<Edge>                EdgeContainer;
  typedef utils::ptr_vector<Triangle>            TriangleContainer;
  typedef utils::ptr_vector<Quad>                QuadContainer;
  typedef std::vector<PtrData>                   DataContainer;
  typedef utils::ptr_vector<PropertyContainer>   PropertyContainerContainer;
  typedef std::vector<std::pair<double, double>> BoundingBox;

  /**
   * @brief Resets the internal geometry ID counter to start from anew.
   *
   * This method is only used for test cases.
   */
  static void resetGeometryIDsGlobally();

  /**
   * @brief Constructor.
   *
   * @param name [IN] Unique name of the mesh.
   * @param flipNormals [IN] Inverts the standard direction of normals.
   */
  Mesh (
    const std::string& name,
    int                dimensions,
    bool               flipNormals );

  /**
   * @brief Destructor, deletes created objects.
   */
  virtual ~Mesh();

  /**
   * @brief Returns group object with all Triangle, Edge, Vertex objects.
   */
  const Group& content();

  /**
   * @brief Returns modifieable container holding all vertices.
   */
  VertexContainer& vertices();

  /**
   * @brief Returns const container holding all vertices.
   */
  const VertexContainer& vertices() const;

  /**
   * @brief Returns modifieable container holding all edges.
   */
  EdgeContainer& edges();

  /**
   * @brief Returns const container holding all edges.
   */
  const EdgeContainer& edges() const;

  /**
   * @brief Returns modifieable container holding all triangles.
   */
  TriangleContainer& triangles();

  /**
   * @brief Returns const container holding all triangles.
   */
  const TriangleContainer& triangles() const;

  /**
   * @brief Returns modifieable container holding all quads.
   */
  QuadContainer& quads();

  /**
   * @brief Returns const container holding all quads.
   */
  const QuadContainer& quads() const;

  PropertyContainerContainer& propertyContainers();

  const PropertyContainerContainer& propertyContainers() const;

  int getDimensions() const;

  /**
   * @brief Registers the listener to be notified when the mesh changes.
   */
  void addListener ( MeshListener& listener );

  template<typename VECTOR_T>
  Vertex& createVertex ( const VECTOR_T& coords )
  {
    assertion2(coords.size() == _dimensions, coords.size(), _dimensions);
    Vertex* newVertex = new Vertex(coords, _manageVertexIDs.getFreeID(), *this);
    newVertex->addParent(*this);
    _content.add(newVertex);
    return *newVertex;
  }

  /**
   * @brief Creates and initializes an Edge object.
   *
   * @param vertexOne [IN] Reference to first Vertex defining the Edge.
   * @param vertexTwo [IN] Reference to second Vertex defining the Edge.
   * @param meshIsParent [IN] If true, the mesh is set as Property parent.
   */
  Edge& createEdge (
    Vertex& vertexOne,
    Vertex& vertexTwo );

  /**
   * @brief Creates and initializes a Triangle object.
   *
   * @param edgeOne [IN] Reference to first edge defining the Triangle.
   * @param edgeTwo [IN] Reference to second edge defining the Triangle.
   * @param edgeThree [IN] Reference to third edge defining the Triangle.
   * @param meshIsParent [IN] If true, the mesh is set as Property parent.
   */
  Triangle& createTriangle (
    Edge& edgeOne,
    Edge& edgeTwo,
    Edge& edgeThree );

  /**
   * @brief Creates and initializes a Quad object.
   *
   * @param edgeOne [IN] Reference to first edge defining the Quad.
   * @param edgeTwo [IN] Reference to second edge defining the Quad.
   * @param edgeThree [IN] Reference to third edge defining the Quad.
   * @param edgeFour [IN] Reference to fourth edge defining the Quad.
   * @param meshIsParent [IN] If true, the mesh is set as Property parent.
   */
  Quad& createQuad (
    Edge& edgeOne,
    Edge& edgeTwo,
    Edge& edgeThree,
    Edge& edgeFour);

  /**
   * @brief Creates and initializes a PropertyContainer object.
   *
   * A PropertyContainer can be used to serve as Property parent for other
   * objects in the mesh. It is used, to define multiple geometry IDs for
   * Vertices and Edges lying at the intersection of regions with different
   * geometry IDs.
   */
  PropertyContainer& createPropertyContainer();

  PtrData& createData (
    const std::string& name,
    int                dimension );

  const DataContainer& data() const;

  const PtrData& data ( int dataID ) const;

  PropertyContainer& getPropertyContainer ( const std::string subIDName );

  const std::string& getName() const;

  bool isFlipNormals() const;

  void setFlipNormals ( bool flipNormals );

  /**
   * @brief Associates a new geometry ID to the mesh.
   *
   * The full name of the ID is given by:
   * getName() + "-" + subIDNamePostfix
   *
   * @return Newly created property container holding the new ID.
   */
  PropertyContainer& setSubID ( const std::string& subIDNamePostfix );

  /**
   * @brief Returns all used geometry IDs paired with their names
   */
  const std::map<std::string,int>& getNameIDPairs();

  /**
   * @brief Returns the geometry ID corresponding to the given name.
   */
  int getID ( const std::string& name ) const;

  /**
   * @brief Returns the base ID of the mesh.
   */
  int getID() const;

  /**
   * @brief Allocates memory for the vertex data values.
   */
  void allocateDataValues();

  /**
   * @brief Necessary before any geom. operations can be performed on the mesh.
   *
   * If no edges (in 2d) or triangles/quads (in 3d) are
   * given, no normals are computed in order to avoid dividing by zero on
   * normalization of the vertex normals.
   *
   * Circumcircles of edges and triangles are computed.
   */
  void computeState();

  /**
   * @brief Removes all mesh elements and data values (does not remove data).
   *
   * A mesh element is a
   * - vertex
   * - edge
   * - triangle
   * - property container
   */
  void clear();

  /**
   * @brief Notifies all MeshListeners that the mesh has changed.
   */
  void notifyListeners();

  std::map<int,std::vector<int> >& getVertexDistribution(){
    return _vertexDistribution;
  }

  int getGlobalNumberOfVertices(){
    return _globalNumberOfVertices;
  }

  void setGlobalNumberOfVertices(int globalNumberOfVertices){
    _globalNumberOfVertices = globalNumberOfVertices;
  }

  void addMesh(Mesh& deltaMesh);

  /**
   * @brief Returns the bounding box of the mesh.
   *
   * BoundingBox is a vector of pairs (min, max), one pair for each dimension.
   * computeState() has to be called after setting the mesh.
   */
  const BoundingBox getBoundingBox() const;

  /**
   * @brief Returns the Center Of Gravity of the mesh
   * 
   * Returns a vector of doubles, size d, each dimension computed as
   * cog =  (max - min) / 2 + min
   */
  const std::vector<double> getCOG() const;
  
private:

  // @brief Logging device.
  static tarch::logging::Log _log;

  // @brief Provides unique IDs for all geometry objects
  static utils::ManageUniqueIDs* _managerPropertyIDs;

  // @brief Name of the mesh.
  std::string _name;

  // @brief Dimension of mesh.
  int _dimensions;

  // @brief Flag for flipping normals direction.
  bool _flipNormals;

  // @brief Holds all mesh names and the corresponding IDs belonging to the mesh.
  std::map<std::string,int> _nameIDPairs;

  // @brief Holds vertices, edges, and triangles.
  Group _content;

  // @brief All property containers created by the mesh.
  PropertyContainerContainer _propertyContainers;

  // @brief Data hold by the vertices of the mesh.
  DataContainer _data;

  utils::ManageUniqueIDs _manageVertexIDs;

  utils::ManageUniqueIDs _manageEdgeIDs;

  utils::ManageUniqueIDs _manageTriangleIDs;

  utils::ManageUniqueIDs _manageQuadIDs;

  // @brief Mesh listeners interested in mesh changes.
  std::list<MeshListener*> _listeners;

  /**
   * @brief Vertex distribution for the master, holding for each slave all vertex IDs.
   */
  std::map<int,std::vector<int> > _vertexDistribution;

  /**
   * @brief Global number of vertices for the master.
   */
  int _globalNumberOfVertices;

  BoundingBox _boundingBox;
};

}} // namespace precice, mesh

#endif /* PRECICE_MESH_MESH_HPP_ */
