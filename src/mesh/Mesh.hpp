#pragma once

#include "mesh/Group.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Data.hpp"
#include "mesh/Vertex.hpp"
#include "utils/PointerVector.hpp"
#include "utils/ManageUniqueIDs.hpp"
#include <map>
#include <list>
#include <vector>
#include <boost/signals2.hpp>

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
class Mesh : public PropertyContainer
{
public:

  using VertexContainer            = utils::ptr_vector<Vertex>;
  using EdgeContainer              = utils::ptr_vector<Edge>;
  using TriangleContainer          = utils::ptr_vector<Triangle>;
  using QuadContainer              = utils::ptr_vector<Quad>;
  using DataContainer              = std::vector<PtrData>;
  using PropertyContainerContainer = utils::ptr_vector<PropertyContainer>;
  using BoundingBox                = std::vector<std::pair<double, double>>;
  using BoundingBoxMap             = std::map<int,BoundingBox>; 

  /// A mapping from rank to used (not necessarily owned) vertex IDs
  using VertexDistribution = std::map<int, std::vector<int>>;

  /// A mapping from remote local ranks to the IDs that must be communicated
  using CommunicationMap = std::map<int, std::vector<int>>;

  /// Signal is emitted when the mesh is changed
  boost::signals2::signal<void(Mesh &)> meshChanged;

  /// Signal is emitted when the mesh is destroyed
  boost::signals2::signal<void(Mesh &)> meshDestroyed;

  /**
   * @brief Resets the internal geometry ID counter to start from anew.
   *
   * This method is only used for test cases.
   */
  static void resetGeometryIDsGlobally();

  /**
   * @brief Constructor.
   *
   * @param[in] name Unique name of the mesh.
   * @param[in] dimensions Dimensionalty of the mesh.
   * @param[in] flipNormals Inverts the standard direction of normals.
   */
  Mesh (
    const std::string& name,
    int                dimensions,
    bool               flipNormals );

  /// Destructor, deletes created objects.
  virtual ~Mesh();

  /// Returns group object with all Triangle, Edge, Vertex objects.
  const Group& content() const;

  /// Returns modifieable container holding all vertices.
  VertexContainer& vertices();

  /// Returns const container holding all vertices.
  const VertexContainer& vertices() const;

  /// Returns modifiable container holding all edges.
  EdgeContainer& edges();

  /// Returns const container holding all edges.
  const EdgeContainer& edges() const;

  /// Returns modifiable container holding all triangles.
  TriangleContainer& triangles();

  /// Returns const container holding all triangles.
  const TriangleContainer& triangles() const;

  /// Returns modifiable container holding all quads.
  QuadContainer& quads();

  /// Returns const container holding all quads.
  const QuadContainer& quads() const;

  PropertyContainerContainer& propertyContainers();

  const PropertyContainerContainer& propertyContainers() const;

  int getDimensions() const;

  template<typename VECTOR_T>
  Vertex& createVertex ( const VECTOR_T& coords )
  {
    PRECICE_ASSERT(coords.size() == _dimensions, coords.size(), _dimensions);
    Vertex* newVertex = new Vertex(coords, _manageVertexIDs.getFreeID());
    newVertex->addParent(*this);
    _content.add(newVertex);
    return *newVertex;
  }

  /**
   * @brief Creates and initializes an Edge object.
   *
   * @param[in] vertexOne Reference to first Vertex defining the Edge.
   * @param[in] vertexTwo Reference to second Vertex defining the Edge.
   */
  Edge& createEdge (
    Vertex& vertexOne,
    Vertex& vertexTwo );

  /**
   * @brief Creates and initializes an Edge object or returns an already existing one.
   *
   * @param[in] vertexOne Reference to first Vertex defining the Edge.
   * @param[in] vertexTwo Reference to second Vertex defining the Edge.
   */
  Edge& createUniqueEdge (
    Vertex& vertexOne,
    Vertex& vertexTwo );

  /**
   * @brief Creates and initializes a Triangle object.
   *
   * @param[in] edgeOne Reference to first edge defining the Triangle.
   * @param[in] edgeTwo Reference to second edge defining the Triangle.
   * @param[in] edgeThree Reference to third edge defining the Triangle.
   */
  Triangle& createTriangle (
    Edge& edgeOne,
    Edge& edgeTwo,
    Edge& edgeThree );

  /**
   * @brief Creates and initializes a Quad object.
   *
   * @param[in] edgeOne Reference to first edge defining the Quad.
   * @param[in] edgeTwo Reference to second edge defining the Quad.
   * @param[in] edgeThree Reference to third edge defining the Quad.
   * @param[in] edgeFour Reference to fourth edge defining the Quad.
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

  PropertyContainer& getPropertyContainer (const std::string & subIDName);

  /// Returns the name of the mesh, as set in the config file.
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

  /// Returns all used geometry IDs paired with their names
  const std::map<std::string,int>& getNameIDPairs();

  /// Returns the geometry ID corresponding to the given name.
  int getID ( const std::string& name ) const;

  /// Returns the base ID of the mesh.
  int getID() const;

  /// Returns true if the given vertexID is valid
  bool isValidVertexID(int vertexID) const;

  /// Returns true if the given edgeID is valid
  bool isValidEdgeID(int edgeID) const;

  /// Allocates memory for the vertex data values.
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

  /// Returns a mapping from rank to used (not necessarily owned) vertex IDs
  VertexDistribution & getVertexDistribution();

  VertexDistribution const & getVertexDistribution() const;

  std::vector<int>& getVertexOffsets();

  const std::vector<int>& getVertexOffsets() const;

  /// Only used for tests
  void setVertexOffsets(std::vector<int> & vertexOffsets);

  int getGlobalNumberOfVertices() const;

  void setGlobalNumberOfVertices(int num);

  /// Returns a vector of connected ranks
  std::vector<int> & getConnectedRanks()
  {
    return _connectedRanks;
  }
  
  /// Returns a mapping from remote local connected ranks to the corresponding vertex IDs
  CommunicationMap & getCommunicationMap()
  {
    return _communicationMap;
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
  
  bool operator==(const Mesh& other) const;
  
  bool operator!=(const Mesh& other) const;

private:

  /// Computes the normals for all primitives.
  void computeNormals();

  /// Computes the boundingBox for the vertices.
  void computeBoundingBox();

  mutable logging::Logger _log{"mesh::Mesh"};

  /// Provides unique IDs for all geometry objects
  static std::unique_ptr<utils::ManageUniqueIDs> _managePropertyIDs;

  /// Name of the mesh.
  std::string _name;

  /// Dimension of mesh.
  int _dimensions;

  /// Flag for flipping normals direction.
  bool _flipNormals;

  /// Holds all mesh names and the corresponding IDs belonging to the mesh.
  std::map<std::string,int> _nameIDPairs;

  /// Holds vertices, edges, and triangles.
  Group _content;

  /// All property containers created by the mesh.
  PropertyContainerContainer _propertyContainers;

  /// Data hold by the vertices of the mesh.
  DataContainer _data;

  utils::ManageUniqueIDs _manageVertexIDs;

  utils::ManageUniqueIDs _manageEdgeIDs;

  utils::ManageUniqueIDs _manageTriangleIDs;

  utils::ManageUniqueIDs _manageQuadIDs;

  /**
   * @brief Vertex distribution for the master, holding for each slave all vertex IDs it owns.
   *
   * For slaves, this data structure is empty and should not be used.
   */
  VertexDistribution _vertexDistribution;

  /// Holds the index of the last vertex for each slave.
  /**
   * The last entry holds the total number of vertices.
   * Needed for the matrix-matrix multiplication of the IMVJ acceleration.
   */
  std::vector<int> _vertexOffsets;


  /**
   * @brief Number of unique vertices for complete distributed mesh.
   *
   * Duplicated vertices are only accounted once.
   */
  int _globalNumberOfVertices = -1;

  /**
   * @brief each rank stores list of connected remote ranks.
   * In the m2n package, this is used to create the initial communication channels.
   */
  std::vector<int> _connectedRanks;

  /**
   * @brief each rank stores list of connected ranks and corresponding vertex IDs here. 
   * In the m2n package, this is used to create the final communication channels.
   */
  CommunicationMap _communicationMap;

  BoundingBox _boundingBox;

};

std::ostream& operator<<(std::ostream& os, const Mesh& q);

}} // namespace precice, mesh
