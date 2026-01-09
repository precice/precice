#pragma once

#include <Eigen/Core>
#include <deque>
#include <iosfwd>
#include <list>
#include <map>
#include <string>
#include <string_view>
#include <vector>

#include "logging/Logger.hpp"
#include "mesh/BoundingBox.hpp"
#include "mesh/Data.hpp"
#include "mesh/Edge.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Tetrahedron.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/Vertex.hpp"
#include "precice/impl/Types.hpp"
#include "query/Index.hpp"
#include "utils/ManageUniqueIDs.hpp"
#include "utils/assertion.hpp"

namespace precice::mesh {

/**
 * @brief Container and creator for meshes.
 *
 * A Mesh can consist of Vertices, Edges, Triangles, and nested Meshes.
 * It provides functionality to conveniently create those objects.
 *
 * In addition to creating the topological information of a mesh, the Mesh class
 * can also be used to add data to the vertices of a mesh, such as function values or also gradient values.
 *
 * Usage example: precice::mesh::tests::MeshTest::testDemonstration()
 */
class Mesh {
public:
  using VertexContainer   = std::deque<Vertex>;
  using EdgeContainer     = std::deque<Edge>;
  using TriangleContainer = std::deque<Triangle>;
  using TetraContainer    = std::deque<Tetrahedron>;
  using DataContainer     = std::vector<PtrData>;
  using BoundingBoxMap    = std::map<int, BoundingBox>;

  /// A mapping from rank to used (not necessarily owned) vertex IDs
  using VertexDistribution = std::map<Rank, std::vector<VertexID>>;

  /// A mapping from remote local ranks to the IDs that must be communicated
  using CommunicationMap = std::map<Rank, std::vector<VertexID>>;
  using ConnectionMap    = CommunicationMap; // until we decide on a name

  using VertexOffsets = std::vector<int>;

  /// Use if the id of the mesh is not necessary
  static constexpr MeshID MESH_ID_UNDEFINED{-1};

  /**
   * @brief Constructor.
   *
   * @param[in] name Unique name of the mesh.
   * @param[in] dimensions Dimensionalty of the mesh.
   * @param[in] id The id of this mesh
   */
  Mesh(
      std::string name,
      int         dimensions,
      MeshID      id,
      bool        isJustInTime = false);

  /// Mutable access to a vertex by VertexID
  Vertex &vertex(VertexID id);

  /// Const access to a vertex by VertexID
  const Vertex &vertex(VertexID id) const;

  /// Returns modifieable container holding all vertices.
  VertexContainer &vertices();

  /// Returns const container holding all vertices.
  const VertexContainer &vertices() const;

  /// Returns the number of vertices
  std::size_t nVertices() const;

  /// Does the mesh contain any vertices?
  bool empty() const
  {
    return _vertices.empty();
  }

  /// Returns modifiable container holding all edges.
  EdgeContainer &edges();

  /// Returns const container holding all edges.
  const EdgeContainer &edges() const;

  bool hasEdges() const
  {
    return !_edges.empty();
  }

  /// Returns modifiable container holding all triangles.
  TriangleContainer &triangles();

  /// Returns const container holding all triangles.
  const TriangleContainer &triangles() const;

  bool hasTriangles() const
  {
    return !_triangles.empty();
  }

  /// Returns modifiable container holding all tetrahedra.
  TetraContainer &tetrahedra();

  /// Returns const container holding all tetrahedra.
  const TetraContainer &tetrahedra() const;

  bool hasTetrahedra() const
  {
    return !_tetrahedra.empty();
  }

  bool hasConnectivity() const
  {
    return hasEdges() || hasTriangles() || hasTetrahedra();
  }

  int getDimensions() const;

  /// Creates and initializes a Vertex object.
  Vertex &createVertex(const Eigen::Ref<const Eigen::VectorXd> &coords);

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

  /**
   * @brief Creates and initializes a Triangle object.
   *
   * @param[in] vertexOne Reference to first edge defining the Triangle.
   * @param[in] vertexTwo Reference to second edge defining the Triangle.
   * @param[in] vertexThree Reference to third edge defining the Triangle.
   */
  Triangle &createTriangle(
      Vertex &vertexOne,
      Vertex &vertexTwo,
      Vertex &vertexThree);

  /**
   * @brief Creates and initializes a Tetrahedron object.
   *
   * @param[in] vertexOne Reference to first vertex defining the Tetrahedron.
   * @param[in] vertexTwo Reference to second vertex defining the Tetrahedron.
   * @param[in] vertexThree Reference to third vertex defining the Tetrahedron.
   * @param[in] vertexFour Reference to fourth vertex defining the Tetrahedron.
   */
  Tetrahedron &createTetrahedron(
      Vertex &vertexOne,
      Vertex &vertexTwo,
      Vertex &vertexThree,
      Vertex &vertexFour);

  /// Create only data for vertex
  PtrData &createData(const std::string                 &name,
                      int                                dimension,
                      DataID                             id,
                      int                                waveformDegree = time::Time::DEFAULT_WAVEFORM_DEGREE,
                      std::vector<std::optional<double>> lowerBound     = std::vector<std::optional<double>>(3),
                      std::vector<std::optional<double>> upperBound     = std::vector<std::optional<double>>(3));

  /// Allows access to all data
  const DataContainer &data() const;

  /// Returns whether Mesh has Data with the matchingID
  bool hasDataID(DataID dataID) const;

  /// Returns the data with the matching ID
  const PtrData &data(DataID dataID) const;

  /// Returns whether Mesh has Data with the dataName
  bool hasDataName(std::string_view dataName) const;

  /// Returns the names of all available data
  std::vector<std::string> availableData() const;

  /// Returns the data with the matching name
  const PtrData &data(std::string_view dataName) const;

  /// Returns the name of the mesh, as set in the config file.
  const std::string &getName() const;

  /// Returns the base ID of the mesh.
  MeshID getID() const;

  /// Returns true if the given vertexID is valid
  bool isValidVertexID(VertexID vertexID) const;

  /// Allocates memory for the vertex data values and corresponding gradient values.
  void allocateDataValues(); //@todo Redesign mapping and remove this function. See https://github.com/precice/precice/issues/1651.

  /// Computes the boundingBox for the vertices.
  void computeBoundingBox();

  /**
   * @brief Removes all mesh elements and data values (does not remove data or the bounding boxes).
   *
   * A mesh element is a
   * - vertex
   * - edge
   * - triangle
   */
  void clear();

  /// Clears the partitioning information
  void clearPartitioning();

  /// Clears all data stamples
  void clearDataStamples();

  void setVertexDistribution(VertexDistribution vd)
  {
    PRECICE_ASSERT(std::all_of(vd.begin(), vd.end(), [](const auto &p) { return std::is_sorted(p.second.begin(), p.second.end()); }));
    _vertexDistribution = std::move(vd);
  }

  /// Returns a mapping from rank to used (not necessarily owned) vertex IDs
  const VertexDistribution &getVertexDistribution() const
  {
    return _vertexDistribution;
  }

  const VertexOffsets &getVertexOffsets() const
  {
    return _vertexOffsets;
  }

  /// checks if the given ranks partition is empty
  bool isPartitionEmpty(Rank rank) const;

  /// Only used for tests
  void setVertexOffsets(VertexOffsets vertexOffsets)
  {
    _vertexOffsets = std::move(vertexOffsets);
  }

  int getGlobalNumberOfVertices() const
  {
    return _globalNumberOfVertices;
  }

  void setGlobalNumberOfVertices(int num)
  {
    _globalNumberOfVertices = num;
  }

  // Get the data of owned vertices for given data ID
  Eigen::VectorXd getOwnedVertexData(const Eigen::VectorXd &values);

  // Tag all the vertices
  void tagAll();

  /// Returns a vector of connected ranks
  const std::vector<Rank> &getConnectedRanks() const
  {
    return _connectedRanks;
  }

  /// Returns a vector of connected ranks
  void setConnectedRanks(std::vector<Rank> ranks)
  {
    _connectedRanks = std::move(ranks);
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
   * Note that the bounding box does not necessarily need to match the bounding
   * box of the contained vertices of the mesh. Examples are the direct mesh access
   * or a computation of the bounding box before applying additional filtering
   * during the repartitioning. Note that mesh::clear doesn't clear the underlying
   * bounding box (which is potentially a user input when using and defining
   * direct mesh access)
   */
  const BoundingBox &getBoundingBox() const;

  void expandBoundingBox(const BoundingBox &bounding_box);

  bool operator==(const Mesh &other) const;

  bool operator!=(const Mesh &other) const;

  /// Call preprocess() before index() to ensure correct projection handling
  const query::Index &index() const
  {
    return _index;
  }

  /// Call preprocess() before index() to ensure correct projection handling
  query::Index &index()
  {
    return _index;
  }

  bool isJustInTime() const
  {
    return _isJustInTime;
  }

  /**
   * Removes all duplicates and generates implicit primitives.
   *
   * This needs to be called for correct projection handling.
   *
   * @see removeDuplicates()
   * @see generateImplictPrimitives()
   */
  void preprocess();

private:
  mutable logging::Logger _log{"mesh::Mesh"};

  /// Name of the mesh.
  std::string _name;

  /// Dimension of mesh.
  int _dimensions;

  /// The ID of this mesh.
  MeshID _id;

  /// Holds vertices, edges, triangles and tetrahedra.
  VertexContainer   _vertices;
  EdgeContainer     _edges;
  TriangleContainer _triangles;
  TetraContainer    _tetrahedra;

  /// Data hold by the vertices of the mesh.
  DataContainer _data;

  /**
   * @brief Vertex distribution for the primary rank, holding for each secondary rank all vertex IDs it owns.
   *
   * For secondary ranks, this data structure is empty and should not be used.
   */
  VertexDistribution _vertexDistribution;

  /// Holds the index of the last vertex for each rank.
  /**
   * The last entry holds the total number of vertices.
   * Needed for the matrix-matrix multiplication of the IMVJ acceleration.
   */
  VertexOffsets _vertexOffsets;

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
  std::vector<Rank> _connectedRanks;

  /**
   * @brief each rank stores list of connected ranks and corresponding vertex IDs here.
   * In the m2n package, this is used to create the final communication channels.
   */
  CommunicationMap _communicationMap;

  /// for just-in-time mapping, we need an artificial mesh, which we can use
  bool _isJustInTime = false;

  BoundingBox _boundingBox;

  query::Index _index;

  /// Removes all duplicate connectivity.
  void removeDuplicates();

  /** Generates implicit primitives for correct projection handling
   *
   * removeDuplicates() should be called first to avoid unnecessary filtering.
   */
  void generateImplictPrimitives();
};

std::ostream &operator<<(std::ostream &os, const Mesh &q);

} // namespace precice::mesh
