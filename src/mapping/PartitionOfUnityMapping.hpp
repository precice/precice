#pragma once

#include <Eigen/Core>
#include <numeric>

#include "com/Communication.hpp"
#include "io/ExportVTU.hpp"
#include "mapping/CreatePartitioning.hpp"
#include "mapping/Partition.hpp"
#include "mapping/RadialBasisFctBaseMapping.hpp"
#include "mesh/Filter.hpp"
#include "precice/types.hpp"
#include "query/Index.hpp"
#include "utils/Event.hpp"
#include "utils/IntraComm.hpp"

namespace precice {
extern bool syncMode;

namespace mapping {

/**
 * @brief Mapping using partition of unity decomposition strategies
 */
template <typename RADIAL_BASIS_FUNCTION_T>
class PartitionOfUnityMapping : public Mapping {
public:
  /**
   * @brief Constructor.
   *
   * @param[in] constraint Specifies mapping to be consistent or conservative.
   * @param[in] dimension Dimensionality of the meshes
   * @param[in] parameter shape parameter or support radius of the interpolation RBF
   * @param[in] verticesPerPartition Estimate for the number of vertices to be clustered together
   * @param[in] relativeOverlap Overlap between partitions, where 1 would correspnd to a complete overlap, 0 to distance of 2 x radius between partitions
   */
  PartitionOfUnityMapping(
      Mapping::Constraint constraint,
      int                 dimension,
      double              parameter,
      Polynomial          polynomial,
      unsigned int        verticesPerPartition,
      double              relativeOverlap,
      bool                projectToInput);

  /// Computes the mapping coefficients from the in- and output mesh.
  virtual void computeMapping() override;

  /// Removes a computed mapping.
  virtual void clear() override;

  /// tag the vertices required for the mapping
  virtual void tagMeshFirstRound() override;

  /// nothing to do here
  virtual void tagMeshSecondRound() override;

private:
  precice::logging::Logger _log{"mapping::PartitionOfUnityMapping"};

  std::vector<Partition<RADIAL_BASIS_FUNCTION_T>> _partitions;

  // Shape parameter or support radius for the RBF interpolant,
  // only required for the Partition instantiation
  const double _parameter;

  // Input parameter
  const unsigned int _verticesPerPartition;
  const double       _relativeOverlap;
  const bool         _projectToInput;

  // Derived parameter
  double averagePartitionRadius = 0;

  /// true if the mapping along some axis should be ignored
  /// has currently only dim x false entries, as integrated polynomials are irrelevant
  std::vector<bool> _deadAxis;

  Polynomial _polynomial;

  // Holds the output vertex -> partition association. Outer vector has the size of the output mesh and inner vector size of the associated partitions
  std::vector<std::vector<VertexID>> _partMap;

  /// @copydoc Mapping::mapConservative
  virtual void mapConservative(DataID inputDataID, DataID outputDataID) override;

  /// @copydoc Mapping::mapConsistent
  virtual void mapConsistent(DataID inputDataID, DataID outputDataID) override;

  void exportPartitionCentersAsVTU(mesh::Mesh &centers);
};

// --------------------------------------------------- HEADER IMPLEMENTATIONS

template <typename RADIAL_BASIS_FUNCTION_T>
PartitionOfUnityMapping<RADIAL_BASIS_FUNCTION_T>::PartitionOfUnityMapping(
    Mapping::Constraint constraint,
    int                 dimension,
    double              parameter,
    Polynomial          polynomial,
    unsigned int        verticesPerPartition,
    double              relativeOverlap,
    bool                projectToInput)
    : Mapping(constraint, dimension),
      _parameter(parameter), _verticesPerPartition(verticesPerPartition), _relativeOverlap(relativeOverlap), _projectToInput(projectToInput), _polynomial(polynomial)
{
  PRECICE_ASSERT(this->getDimensions() <= 3);
  PRECICE_CHECK(_relativeOverlap < 1, "The relative overlap has to be smaller than one.");
  PRECICE_CHECK(_verticesPerPartition > 0, "The number of vertices per partition has to be greater zero.");

  if (isScaledConsistent()) {
    setInputRequirement(Mapping::MeshRequirement::FULL);
    setOutputRequirement(Mapping::MeshRequirement::FULL);
  } else {
    setInputRequirement(Mapping::MeshRequirement::VERTEX);
    setOutputRequirement(Mapping::MeshRequirement::VERTEX);
  }
  // Simply set each axis to be active according to the space dimension.
  // We could even set each axis to be true, which would just lead to 'zero' columns in the polynomial QR.
  // Integrated polynomials are anyway not supported with PoU
  _deadAxis.clear();
  std::fill_n(std::back_inserter(_deadAxis), getDimensions(), false);
}

template <typename RADIAL_BASIS_FUNCTION_T>
void PartitionOfUnityMapping<RADIAL_BASIS_FUNCTION_T>::tagMeshFirstRound()
{
  PRECICE_TRACE();
  mesh::PtrMesh filterMesh, outMesh;
  if (this->hasConstraint(Mapping::CONSERVATIVE)) {
    filterMesh = this->output(); // remote
    outMesh    = this->input();  // local
  } else {
    filterMesh = this->input();  // remote
    outMesh    = this->output(); // local
  }

  if (outMesh->vertices().empty())
    return; // Ranks not at the interface should never hold interface vertices

  // TODO: Check again the tagging in combination with the partition construction (which mesh to use)
  // In order to construct the local partitions, we need all vertices with a distance of 2 x radius,
  // as the relevant partitions centers have a maximum distance of radius, and the proper construction of the
  // interpolant requires all vertices with a distance of radius from the center.
  auto bb = outMesh->getBoundingBox();

  if (averagePartitionRadius == 0)
    averagePartitionRadius = partitioner::estimatePartitionRadius(_verticesPerPartition, filterMesh, bb);

  // @TODO: This assert is not completely right, as it checks all dimensions for non-emptyness (which might not be the case).
  // However, with the current BB implementation, the expandBy function will just do nothing.
  PRECICE_ASSERT(!bb.empty());
  PRECICE_ASSERT(averagePartitionRadius > 0);
  // Now we extend the bounding box by the radius
  bb.expandBy(1 * averagePartitionRadius);

  // ... and tag all affected vertices
  auto verticesNew = filterMesh->index().getVerticesInsideBox(bb);

  std::for_each(verticesNew.begin(), verticesNew.end(), [&filterMesh](VertexID v) { filterMesh->vertices()[v].tag(); });
}

template <typename RADIAL_BASIS_FUNCTION_T>
void PartitionOfUnityMapping<RADIAL_BASIS_FUNCTION_T>::tagMeshSecondRound()
{
  // Probably nothing to be done here. There is no global ownership for matrix entries required and we tag all potentially locally relevant vertices already in the first round.
}

template <typename RADIAL_BASIS_FUNCTION_T>
void PartitionOfUnityMapping<RADIAL_BASIS_FUNCTION_T>::computeMapping()
{
  PRECICE_TRACE();

  precice::utils::Event e("map.pou.computeMapping.From" + this->input()->getName() + "To" + this->output()->getName(), precice::syncMode);

  // Recompute the whole partitioning
  this->clear();

  mesh::PtrMesh inMesh;
  mesh::PtrMesh outMesh;
  if (this->hasConstraint(Mapping::CONSERVATIVE)) {
    inMesh  = this->output();
    outMesh = this->input();
  } else { // Consistent or scaled consistent
    inMesh  = this->input();
    outMesh = this->output();
  }

  // Step 1: get a tentative partitioning consisting of centers and a radius from one of the available algorithms
  auto [averagePartitionRadius_, centerCandidates] = partitioner::createUniformBlockPartitioning(inMesh, outMesh, _relativeOverlap, _verticesPerPartition, _projectToInput);

  averagePartitionRadius = averagePartitionRadius_;
  PRECICE_ASSERT(averagePartitionRadius > 0);

  // Step 2: check, which of the resulting partitions are be non-empty (in term sof the output vertices) and register the partition centers in a mesh
  mesh::Mesh centerMesh("partitionCentersMesh", this->getDimensions(), mesh::Mesh::MESH_ID_UNDEFINED);
  auto &     meshVertices = centerMesh.vertices();

  for (const auto &c : centerCandidates) {
    // We cannot simply take the vertex from the container, as the ID needs to match the partition ID
    // That's required for the indexing and asserted below
    mesh::Vertex                       center(c.getCoords(), meshVertices.size());
    Partition<RADIAL_BASIS_FUNCTION_T> p(inMesh->getDimensions(), center, averagePartitionRadius, _parameter, _deadAxis, _polynomial, _verticesPerPartition, inMesh, outMesh);

    // Consider only non-empty partitions
    if (!p.isEmpty()) {
      PRECICE_ASSERT(center.getID() == _partitions.size(), center.getID(), _partitions.size());
      meshVertices.emplace_back(std::move(center));
      _partitions.emplace_back(std::move(p));
    }
  }
  // Log the average number of resulting partitions
  unsigned int nVertices   = std::accumulate(_partitions.begin(), _partitions.end(), static_cast<unsigned int>(0), [](auto &acc, auto &val) { return acc += val.getNumberOfInputVertices(); });
  unsigned int maxVertices = std::max_element(_partitions.begin(), _partitions.end(), [](auto &v1, auto &v2) { return v1.getNumberOfInputVertices() < v2.getNumberOfInputVertices(); })->getNumberOfInputVertices();
  unsigned int minVertices = std::min_element(_partitions.begin(), _partitions.end(), [](auto &v1, auto &v2) { return v1.getNumberOfInputVertices() < v2.getNumberOfInputVertices(); })->getNumberOfInputVertices();
  PRECICE_INFO("Average number of vertices per partition {}", nVertices / _partitions.size());
  PRECICE_INFO("Maximum number of vertices per partition {}", maxVertices);
  PRECICE_INFO("Minimum number of vertices per partition {}", minVertices);
  PRECICE_INFO("Number of total partitions (final): {}", _partitions.size());

  // Log a bounding box of the center mesh
  centerMesh.computeBoundingBox();
  PRECICE_INFO("Bounding Box of partition centers {}", centerMesh.getBoundingBox());

  // Step 3: index the partitions / the center mesh in order to create the output vertex -> partition association
  query::Index partitionIndex(centerMesh);
  // Find all partitions the output vertex lies in, i.e., find all centers which have the distance of a partition radius
  // Here: VertexID = partitionID
  // This could also be done on-the-fly in the map data phase, which would require to make the mesh as well as the indexTree member variables.
  PRECICE_ASSERT(_partMap.empty());
  PRECICE_DEBUG("Computing partition-vertex association");
  for (const auto &v : outMesh->vertices()) {
    auto indices = partitionIndex.getVerticesInsideBox(v, averagePartitionRadius);
    if (indices.size() == 0) {
      PRECICE_WARN("Output vertex {} could not be assigned to a partition. This means that the meshes probably do not match well geometry-wise.", v.getCoords());
      // TODO: Think about a proper way to handle this case, maybe set all radii to distance(v, closestvertex)?
      indices.emplace_back(partitionIndex.getClosestVertex(v.getCoords()).index);
    }
    _partMap.emplace_back(indices);
  }

  // Consistency check
  PRECICE_ASSERT(_partMap.size() == outMesh->vertices().size());

// Add a VTK export for visualization purposes
#ifndef NDEBUG
  exportPartitionCentersAsVTU(centerMesh);
#endif

  // Compute the weigths
  precice::utils::Event e_exec("map.pou.computeMapping.computeWeights", precice::syncMode);

  // Iterate over all vertices and update the output data
  for (VertexID v = 0; v < outMesh->vertices().size(); ++v) {
    // 1. Find all partitions the output vertex lies in, i.e., find all centers which have the distance of a partition radius
    std::vector<VertexID> partitionIDs = _partMap[v];
    const auto &          vertex       = outMesh->vertices()[v];

    PRECICE_ASSERT(vertex.getID() == v);
    PRECICE_ASSERT(partitionIDs.size() > 0, "No partition found for vertex v ", v);

    // 2. In each partition, gather the weights
    std::vector<double> weights(partitionIDs.size());
    std::transform(partitionIDs.cbegin(), partitionIDs.cend(), weights.begin(), [&](const auto &ids) { return _partitions[ids].computeWeight(vertex); });
    double weightSum = std::accumulate(weights.begin(), weights.end(), static_cast<double>(0.));
    // TODO: This covers the edge case of vertices being at the edge of (several) partitions
    // In case the sum is equal to zero, we assign equal weights for all partitions
    if (!(weightSum > 0)) {
      PRECICE_ASSERT(weights.size() > 0);
      std::for_each(weights.begin(), weights.end(), [&weights](auto &w) { w = 1 / weights.size(); });
      weightSum = 1;
    }
    PRECICE_DEBUG("Weight sum {}", weightSum);
    PRECICE_DEBUG("Partitions {}", partitionIDs);
    PRECICE_DEBUG("V coords {}", vertex.getCoords());
    PRECICE_ASSERT(weightSum > 0);

    for (unsigned int i = 0; i < partitionIDs.size(); ++i) {
      _partitions[partitionIDs[i]].setNormalizedWeight(weights[i] / weightSum, vertex.getID());
    }
  }

  // Set the computedMapping flag
  this->_hasComputedMapping = true;
  PRECICE_DEBUG("Compute Mapping is Completed.");
}

template <typename RADIAL_BASIS_FUNCTION_T>
void PartitionOfUnityMapping<RADIAL_BASIS_FUNCTION_T>::clear()
{
  PRECICE_TRACE();
  _partitions.clear();
  _partMap.clear();
  // TODO: Don't reset this here
  averagePartitionRadius    = 0;
  this->_hasComputedMapping = false;
}

template <typename RADIAL_BASIS_FUNCTION_T>
void PartitionOfUnityMapping<RADIAL_BASIS_FUNCTION_T>::mapConservative(DataID inputDataID, DataID outputDataID)
{
  precice::utils::Event e("map.pou.mapData.From" + input()->getName() + "To" + output()->getName(), precice::syncMode);

  PRECICE_TRACE(inputDataID, outputDataID);
  // Execute the actual mapping evaluation in all partitions

  // 1. Reset the output data values as we need to accumulate data across partitions later on
  output()->data(outputDataID)->values().setZero();

  PRECICE_ASSERT(_partMap.size() == this->input()->vertices().size(), _partMap.size(), this->input()->vertices().size());

  // 2. Iterate over all partitoins and accumulate the result
  std::for_each(_partitions.begin(), _partitions.end(), [&](auto &p) { p.mapConservative(input()->data(inputDataID),
                                                                                         output()->data(outputDataID)); });
  // Set mapping finished
  std::for_each(_partitions.begin(), _partitions.end(), [&](auto &p) { p.setMappingFinished(); });
}

template <typename RADIAL_BASIS_FUNCTION_T>
void PartitionOfUnityMapping<RADIAL_BASIS_FUNCTION_T>::mapConsistent(DataID inputDataID, DataID outputDataID)
{
  precice::utils::Event e("map.pou.mapData.From" + input()->getName() + "To" + output()->getName(), precice::syncMode);

  // More detailed measurements
  PRECICE_TRACE(inputDataID, outputDataID);

  // 1. Reset the output data values as we need to accumulate data across partitions later on
  output()->data(outputDataID)->values().setZero();

  // 2. Execute the actual mapping evaluation in all partitions and acccumulate the data
  std::for_each(_partitions.begin(), _partitions.end(), [&](auto &p) { p.mapConsistent(input()->data(inputDataID),
                                                                                       output()->data(outputDataID)); });

  // Set mapping finished
  std::for_each(_partitions.begin(), _partitions.end(), [&](auto &p) { p.setMappingFinished(); });
}

template <typename RADIAL_BASIS_FUNCTION_T>
void PartitionOfUnityMapping<RADIAL_BASIS_FUNCTION_T>::exportPartitionCentersAsVTU(mesh::Mesh &centerMesh)
{
  auto dataRadius      = centerMesh.createData("radius", 1, -1);
  auto dataCardinality = centerMesh.createData("partition-size", 1, -1);
  centerMesh.allocateDataValues();
  dataRadius->values().fill(averagePartitionRadius);
  for (unsigned int i = 0; i < _partitions.size(); ++i) {
    dataCardinality->values()[i] = static_cast<double>(_partitions[i].getNumberOfInputVertices());
  }

  // We have to create the global offsets in order to export things in parallel
  if (utils::IntraComm::isSecondary()) {
    // send number of vertices
    PRECICE_DEBUG("Send number of vertices: {}", centerMesh.vertices().size());
    int numberOfVertices = centerMesh.vertices().size();
    utils::IntraComm::getCommunication()->send(numberOfVertices, 0);

    // receive vertex offsets
    mesh::Mesh::VertexOffsets vertexOffsets;
    utils::IntraComm::getCommunication()->broadcast(vertexOffsets, 0);
    PRECICE_DEBUG("My vertex offsets: {}", vertexOffsets);
    PRECICE_ASSERT(centerMesh.getVertexOffsets().empty());
    centerMesh.setVertexOffsets(std::move(vertexOffsets));
  } else if (utils::IntraComm::isPrimary()) {

    mesh::Mesh::VertexOffsets vertexOffsets(utils::IntraComm::getSize());
    vertexOffsets[0] = centerMesh.vertices().size();

    // receive number of secondary vertices and fill vertex offsets
    for (int secondaryRank : utils::IntraComm::allSecondaryRanks()) {
      int numberOfSecondaryRankVertices = -1;
      utils::IntraComm::getCommunication()->receive(numberOfSecondaryRankVertices, secondaryRank);
      PRECICE_ASSERT(numberOfSecondaryRankVertices >= 0);
      vertexOffsets[secondaryRank] = numberOfSecondaryRankVertices + vertexOffsets[secondaryRank - 1];
    }

    // broadcast vertex offsets
    PRECICE_DEBUG("My vertex offsets: {}", centerMesh.getVertexOffsets());
    utils::IntraComm::getCommunication()->broadcast(vertexOffsets);
    centerMesh.setVertexOffsets(std::move(vertexOffsets));
  }

  io::ExportVTU exporter;
  exporter.doExport("pouCenters", "exports", centerMesh);
}
} // namespace mapping
} // namespace precice
