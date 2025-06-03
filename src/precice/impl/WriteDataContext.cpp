#include "WriteDataContext.hpp"

namespace precice::impl {

logging::Logger WriteDataContext::_log{"impl::WriteDataContext"};

WriteDataContext::WriteDataContext(mesh::PtrData data,
                                   mesh::PtrMesh mesh)
    : DataContext(data, mesh),
      _writeDataBuffer(data->getDimensions())
{
}

void WriteDataContext::resetBufferedData()
{
  invalidateMappingCacheAndResetData();
  _writeDataBuffer.values.setZero();
  _writeDataBuffer.gradients.setZero();
}

void WriteDataContext::trimAfter(double time)
{
  _providedData->timeStepsStorage().trimAfter(time);

  // reset all toData
  PRECICE_ASSERT(!hasReadMapping(), "Read mapping is not allowed for WriteDataContext.");
  if (hasWriteMapping()) {
    std::for_each(_mappingContexts.begin(), _mappingContexts.end(), [time](auto &context) { context.toData->timeStepsStorage().trimAfter(time); });
  }
}

void WriteDataContext::completeJustInTimeMapping()
{
  PRECICE_TRACE();
  if (justInTimeMapping) {
    // finalize  mapping the data stored in the cache and transfer it to the _writeDataBuffer
    Eigen::Map<Eigen::MatrixXd> map(_writeDataBuffer.values.data(), _providedData->getDimensions(), _writeDataBuffer.values.size() / _providedData->getDimensions());
    justInTimeMapping->completeJustInTimeMapping(*mappingCache, map);
  }
}

void WriteDataContext::writeAndMapValues(::precice::span<const double> coordinates, ::precice::span<const double> values)
{
  PRECICE_TRACE();
  PRECICE_ASSERT(mappingCache);
  PRECICE_CHECK(justInTimeMapping,
                "This participant attempted to write data to mesh \"{}\" using a just-in-time mapping, "
                "but there is no write mapping configured for that mesh. "
                "Perhaps you forgot to define a <mapping:... direction=\"write\" /> or want to use direct access with the API function \"writeData({0},...)\".",
                getMeshName(), getDataName());
  PRECICE_ASSERT((coordinates.size() / getSpatialDimensions()) * getDataDimensions() == values.size());
  PRECICE_ASSERT(_writeDataBuffer.values.data());

  // We forward both, the _writeDataBuffer and the cache to the justInTimeMapping
  Eigen::Map<const Eigen::MatrixXd> coords(coordinates.data(), getSpatialDimensions(), coordinates.size() / getSpatialDimensions());
  Eigen::Map<const Eigen::MatrixXd> inputData(values.data(), getDataDimensions(), coordinates.size() / getDataDimensions());
  Eigen::Map<Eigen::MatrixXd>       localData(_writeDataBuffer.values.data(), getDataDimensions(), getMeshVertexCount());

  // Function to fill the localData
  justInTimeMapping->mapConservativeAt(coords, inputData, *mappingCache, localData);
}

void WriteDataContext::writeValuesIntoDataBuffer(::precice::span<const VertexID> vertices, ::precice::span<const double> values)
{
  PRECICE_ASSERT(vertices.size() * getDataDimensions() == values.size());
  PRECICE_ASSERT(_writeDataBuffer.values.data());

  Eigen::Map<const Eigen::MatrixXd> inputData(values.data(), getDataDimensions(), vertices.size());
  Eigen::Map<Eigen::MatrixXd>       localData(_writeDataBuffer.values.data(), getDataDimensions(), getMeshVertexCount());

  for (int i = 0; i < static_cast<int>(vertices.size()); ++i) {
    PRECICE_ASSERT(vertices[i] < localData.cols());
    localData.col(vertices[i]) = inputData.col(i);
  }
}

void WriteDataContext::writeGradientsIntoDataBuffer(::precice::span<const VertexID> vertices, ::precice::span<const double> gradients)
{
  const auto gradientComponents = getSpatialDimensions() * getDataDimensions();

  PRECICE_ASSERT(gradientComponents * vertices.size() == gradients.size());
  PRECICE_ASSERT(_writeDataBuffer.gradients.data());

  Eigen::Map<const Eigen::MatrixXd> inputGradients(gradients.data(), gradientComponents, vertices.size());
  Eigen::Map<Eigen::MatrixXd>       localGradients(_writeDataBuffer.gradients.data(), gradientComponents, getMeshVertexCount());

  for (int i = 0; i < static_cast<int>(vertices.size()); ++i) {
    PRECICE_ASSERT(vertices[i] < localGradients.cols());
    localGradients.col(vertices[i]) = inputGradients.col(i);
  }
}

void WriteDataContext::resizeBufferTo(int nVertices)
{
  using SizeType = std::remove_cv<decltype(nVertices)>::type;

  // Allocate data values
  const SizeType expectedSize = nVertices * getDataDimensions();
  const auto     actualSize   = static_cast<SizeType>(_writeDataBuffer.values.size());
  const auto     change       = expectedSize - actualSize;

  _writeDataBuffer.values.conservativeResize(expectedSize);
  // Enlarge Buffer
  if (change > 0) {
    _writeDataBuffer.values.tail(change).setZero();
  }
  PRECICE_DEBUG("Data {} now has {} values", getDataName(), _writeDataBuffer.values.size());

  if (!_providedData->hasGradient()) {
    return;
  }

  // Allocate gradient data values
  _writeDataBuffer.gradients.conservativeResize(getSpatialDimensions(), expectedSize);

  // Enlarge Buffer
  if (change > 0) {
    _writeDataBuffer.gradients.rightCols(change).setZero();
  }
  PRECICE_DEBUG("Gradient Data {} now has {} x {} values", getDataName(), _writeDataBuffer.gradients.rows(), _writeDataBuffer.gradients.cols());
}

void WriteDataContext::storeBufferedData(double currentTime)
{
  _providedData->setSampleAtTime(currentTime, _writeDataBuffer);
}

void WriteDataContext::appendMappingConfiguration(MappingContext &mappingContext, const MeshContext &meshContext)
{
  PRECICE_ASSERT(meshContext.mesh->hasDataName(getDataName()));
  mesh::PtrData data = meshContext.mesh->data(getDataName());
  PRECICE_ASSERT(data != _providedData, "Data the write mapping is mapping to needs to be different from _providedData");
  mappingContext.fromData = _providedData;
  mappingContext.toData   = data;
  appendMapping(mappingContext);
  PRECICE_ASSERT(hasWriteMapping());
}

} // namespace precice::impl
