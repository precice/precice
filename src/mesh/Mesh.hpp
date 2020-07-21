#pragma once

#include <Eigen/Core>
#include <boost/signals2.hpp>
#include <deque>
#include <iosfwd>
#include <list>
#include <map>
#include <string>
#include <vector>
#include "logging/Logger.hpp"
#include "mesh/BoundingBox.hpp"
#include "mesh/Data.hpp"
#include "mesh/Edge.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/Vertex.hpp"
#include "utils/ManageUniqueIDs.hpp"
#include "utils/PointerVector.hpp"
#include "utils/assertion.hpp"

namespace precice {
namespace mesh {

/**
 * @brief Container and creator for meshes.
 *
 * A Mesh can consist of Vertices, Edges, Triangles, and nested Meshes.
 * It provides functionality to conveniently create those objects.
 *
 * In addition to creating the topological information of a mesh, the Mesh class
 * can also be used to add data to the vertices of a mesh.
 *
 * Usage example: precice::mesh::tests::MeshTest::testDemonstration()
 */
class Mesh {
public:
  using VertexContainer   = std::deque<Vertex>;
  using EdgeContainer     = std::deque<Edge>;
  using TriangleContainer = std::deque<Triangle>;
  using DataContainer     = std::vector<PtrData>;
  using BoundingBoxMap    = std::map<int, BoundingBox>;

  /// A mapping from rank to used (not necessarily owned) vertex IDs
  using VertexDistribution = std::map<int, std::vector<int>>;

  /// A mapping from remote local ranks to the IDs that must be communicated
  using CommunicationMap = std::map<int, std::vector<int>>;

  /// Signal is emitted when the mesh is changed
  boost::signals2::signal<void(Mesh &)> meshChanged;

  /// Signal is emitted when the mesh is destroyed
  boost::signals2::signal<void(Mesh &)> meshDestroyed;

  /// Use if the id of the mesh is not necessary
  static constexpr int MESH_ID_UNDEFINED{-1};

  /**
   * @brief Constructor.
   *
   * @param[in] name Unique name of the mesh.
   * @param[in] dimensions Dimensionalty of the mesh.
   * @param[in] flipNormals Inverts the standard direction of normals.
   * @param[in] id The id of this mesh
   */
  Mesh(
      const std::string &name,
      int                dimensions,
      bool               flipNormals,
      int                id);

  /// Destructor, deletes created objects.
  ~Mesh();

  /// Returns modifieable container holding all vertices.
  VertexContainer &vertices();

  /// Returns const container holding all vertices.
  const VertexContainer &vertices() const;

  /// Returns modifiable container holding all edges.
  EdgeContainer &edges();

  /// Returns const container holding all edges.
  const EdgeContainer &edges() const;

  /// Returns modifiable container holding all triangles.
  TriangleContainer &triangles();

  /// Returns const container holding all triangles.
  const TriangleContainer &triangles() const;

  int getDimensions() const;

  template <typename VECTOR_T>
  Vertex &createVertex(const VECTOR_T &coords)
  {
    PRECICE_ASSERT(coords.size() == _dimensions, coords.size(), _dimensions);
    _vertices.emplace_back(coords, _manageVertexIDs.getFreeID());
    return _vertices.back();
  }

  /**
   * @brief Creates and initializes an Edge object.
   *
   * @param[in] vertexOne Reference to first Vertex defining the Edge.
   * @param[in] vertexTwo Reference to second Vertex defining the Edge.
   */
  Edge &createEdge(
      Vertex &vertexOne,
      Vertex &vertexTwo);

  /**
   * @brief Creates and initializes an Edge object or returns an already existing one.
   *
   * @param[in] vertexOne Reference to first Vertex defining the Edge.
   * @param[in] vertexTwo Reference to second Vertex defining the Edge.
   */
  Edge &createUniqueEdge(
      Vertex &vertexOne,
      Vertex &vertexTwo);

  /**
   * @brief Creates and initializes a Triangle object.
   *
   * @param[in] edgeOne Reference to first edge defining the Triangle.
   * @param[in] edgeTwo Reference to second edge defining the Triangle.
   * @param[in] edgeThree Reference to third edge defining the Triangle.
   */
  Triangle &createTriangle(
      Edge &edgeOne,
      Edge &edgeTwo,
      Edge &edgeThree);

  PtrData &createData(
      const std::string &name,
      int                dimension);

  const DataContainer &data() const;

  const PtrData &data(int dataID) const;

  /// Returns the name of the mesh, as set in the config file.
  const std::string &getName() const;

  bool isFlipNormals() const;

  void setFlipNormals(bool flipNormals);

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
   * If no edges (in 2d) or triangles(in 3d) are
   * given, no normals are computed in order to avoid dividing by zero on
   * normalization of the vertex normals.
   *
   * Circumcircles of edges and triangles are computed.
   */
  void computeState();

  /// Computes the boundingBox for the vertices.
  void computeBoundingBox();

  /**
   * @brief Removes all mesh elements and data values (does not remove data).
   *
   * A mesh element is a
   * - vertex
   * - edge
   * - triangle
   */
  void clear();

  /// Returns a mapping from rank to used (not necessarily owned) vertex IDs
  VertexDistribution &getVertexDistribution();

  VertexDistribution const &getVertexDistribution() const;

  std::vector<int> &getVertexOffsets();

  const std::vector<int> &getVertexOffsets() const;

  /// Only used for tests
  void setVertexOffsets(std::vector<int> &vertexOffsets);

  int getGlobalNumberOfVertices() const;

  void setGlobalNumberOfVertices(int num);

  // Get the data of owned vertices for given data ID
  Eigen::VectorXd getOwnedVertexData(int dataID);

  // Tag all the vertices
  void tagAll();

  /// Returns a vector of connected ranks
  std::vector<int> &getConnectedRanks()
  {
    return _connectedRanks;
  }

  /// Returns a mapping from remote local connected ranks to the corresponding vertex IDs
  CommunicationMap &getCommunicationMap()
  {
    return _communicationMap;
  }

  void addMesh(Mesh &deltaMesh);

  /**
   * @brief Returns the bounding box of the mesh.
   *
   * BoundingBox is a vector of pairs (min, max), one pair for each dimension.
   * computeState() has to be called after setting the mesh.
   */
  const BoundingBox &getBoundingBox() const;

  bool operator==(const Mesh &other) const;

  bool operator!=(const Mesh &other) const;

private:
  mutable logging::Logger _log{"mesh::Mesh"};

  /// Name of the mesh.
  std::string _name;

  /// Dimension of mesh.
  int _dimensions;

  /// Flag for flipping normals direction.
  bool _flipNormals;

  /// The ID of this mesh.
  int _id;

  /// Holds vertices, edges, and triangles.
  VertexContainer   _vertices;
  EdgeContainer     _edges;
  TriangleContainer _triangles;

  /// Data hold by the vertices of the mesh.
  DataContainer _data;

  utils::ManageUniqueIDs _manageVertexIDs;

  utils::ManageUniqueIDs _manageEdgeIDs;

  utils::ManageUniqueIDs _manageTriangleIDs;

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

std::ostream &operator<<(std::ostream &os, const Mesh &q);

} // namespace mesh
} // namespace precice
