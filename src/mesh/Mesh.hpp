#pragma once

#include "mesh/Group.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Data.hpp"
#include "mesh/Vertex.hpp"
#include "utils/PointerVector.hpp"
#include "utils/ManageUniqueIDs.hpp"
#include <boost/noncopyable.hpp>
#include <map>
#include <list>
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
   * @param[in] name Unique name of the mesh.
   * @param[in] flipNormals Inverts the standard direction of normals.
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
    assertion(coords.size() == _dimensions, coords.size(), _dimensions);
    Vertex* newVertex = new Vertex(coords, _manageVertexIDs.getFreeID(), *this);
    newVertex->addParent(*this);
    _content.add(newVertex);
    return *newVertex;
  }

  /**
   * @brief Creates and initializes an Edge object.
   *
   * @param[in] vertexOne Reference to first Vertex defining the Edge.
   * @param[in] vertexTwo Reference to second Vertex defining the Edge.
   * @param[in] meshIsParent If true, the mesh is set as Property parent.
   */
  Edge& createEdge (
    Vertex& vertexOne,
    Vertex& vertexTwo );

  /**
   * @brief Creates and initializes a Triangle object.
   *
   * @param[in] edgeOne Reference to first edge defining the Triangle.
   * @param[in] edgeTwo Reference to second edge defining the Triangle.
   * @param[in] edgeThree Reference to third edge defining the Triangle.
   * @param[in] meshIsParent If true, the mesh is set as Property parent.
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
   * @param[in] meshIsParent If true, the mesh is set as Property parent.
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
   * @brief collect/compute global distribution information (for parallel runs) that
   * are needed at a local level, e.g. global indices
   */
  void computeDistribution();

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

  std::vector<int>& getVertexOffsets(){
    return _vertexOffsets;
  }

  //only used for tests
  void setVertexOffsets(std::vector<int> vertexOffsets){
    _vertexOffsets = vertexOffsets;
  }

  int getGlobalNumberOfVertices(){
    return _globalNumberOfVertices;
  }

  void setGlobalNumberOfVertices(int num){
    _globalNumberOfVertices = num;
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

  /// Logging device.
  static logging::Logger _log;

  /// Provides unique IDs for all geometry objects
  static std::unique_ptr<utils::ManageUniqueIDs> _managerPropertyIDs;

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

  /// Mesh listeners interested in mesh changes.
  std::list<MeshListener*> _listeners;

  /**
   * @brief Vertex distribution for the master, holding for each slave all vertex IDs it owns.
   * For slaves, this data structure is empty and should not be used.
   */
  std::map<int,std::vector<int> > _vertexDistribution;

  /// Holds the index of the last vertex for each slave.
  /**
   * The last entry holds the total number of vertices.
   * Needed for the matrix-matrix multiplication of the IMVJ postprocessing.
   */
  std::vector<int> _vertexOffsets;


  /**
   * @brief Number of unique vertices for complete distributed mesh.
   * Duplicated vertices are only accounted once.
   */
  int _globalNumberOfVertices;

  BoundingBox _boundingBox;

  /// Sets the globalIndices on all vertices in the mesh
  void setGlobalIndices(const std::vector<int> &globalIndices);

  /// Sets the isOwner property of a vertex i, if ownerVec[i] == 1
  void setOwnerInformation(const std::vector<int> &ownerVec);
};

}} // namespace precice, mesh
