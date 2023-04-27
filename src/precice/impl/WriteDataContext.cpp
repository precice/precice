#include "WriteDataContext.hpp"

#include "utils/EigenHelperFunctions.hpp"

namespace precice::impl {

logging::Logger WriteDataContext::_log{"impl::WriteDataContext"};

WriteDataContext::WriteDataContext(
    mesh::PtrData data,
    mesh::PtrMesh mesh)
    : DataContext(data, mesh)
{
  _writeDataBuffer = time::Sample{Eigen::VectorXd(), Eigen::MatrixXd()};
}

mesh::PtrData WriteDataContext::providedData()
{
  PRECICE_ASSERT(_providedData);
  return _providedData;
}

time::Sample WriteDataContext::writeDataBuffer()
{
  return _writeDataBuffer;
}

void WriteDataContext::writeValuesIntoDataBuffer(int index, double value)
{
  PRECICE_DEBUG("Store value {} at id {}", value, index);
  PRECICE_ASSERT(_writeDataBuffer.values.size() > index, _writeDataBuffer.values.size(), index);
  _writeDataBuffer.values[index] = value;
}

void WriteDataContext::writeGradientIntoDataBuffer(std::vector<int> indices, const Eigen::Map<const Eigen::MatrixXd> gradients)
{
  const auto vertexCount = getDataSize() / getDataDimensions();
  const int  stride      = getDataDimensions();
  PRECICE_ASSERT(providedData() != nullptr);

  for (auto i = 0; i < indices.size(); i++) {
    const auto valueIndex = indices[i];
    PRECICE_CHECK(0 <= valueIndex && valueIndex < vertexCount,
                  "Cannot write gradient data \"{}\" to invalid Vertex ID ({}). Please make sure you only use the results from calls to setMeshVertex/Vertices().",
                  getDataName(), valueIndex);
    _writeDataBuffer.gradient.block(0, stride * valueIndex, getSpatialDimensions(), getDataDimensions()) = gradients.block(0, stride * i, getSpatialDimensions(), getDataDimensions());
  }
}

void WriteDataContext::resizeBufferTo(int nVertices)
{
  using SizeType = std::remove_cv<decltype(nVertices)>::type;

  // Allocate data values
  const SizeType expectedSize = nVertices * getDataDimensions();
  const auto     actualSize   = static_cast<SizeType>(_writeDataBuffer.values.size());
  // Shrink Buffer
  if (expectedSize < actualSize) {
    _writeDataBuffer.values.resize(expectedSize);
  }
  // Enlarge Buffer
  if (expectedSize > actualSize) {
    const auto leftToAllocate = expectedSize - actualSize;
    utils::append(_writeDataBuffer.values, Eigen::VectorXd(Eigen::VectorXd::Zero(leftToAllocate)));
  }
  PRECICE_DEBUG("Data {} now has {} values", getDataName(), _writeDataBuffer.values.size());

  // Allocate gradient data values
  if (_providedData->hasGradient()) {
    const SizeType spaceDimensions = getSpatialDimensions();

    const SizeType expectedColumnSize = expectedSize * getDataDimensions();
    const auto     actualColumnSize   = static_cast<SizeType>(_writeDataBuffer.gradient.cols());

    // Shrink Buffer
    if (expectedColumnSize < actualColumnSize) {
      _writeDataBuffer.gradient.resize(spaceDimensions, expectedColumnSize);
    }

    // Enlarge Buffer
    if (expectedColumnSize > actualColumnSize) {
      const auto columnLeftToAllocate = expectedColumnSize - actualColumnSize;
      utils::append(_writeDataBuffer.gradient, Eigen::MatrixXd(Eigen::MatrixXd::Zero(spaceDimensions, columnLeftToAllocate)));
    }
    PRECICE_DEBUG("Gradient Data {} now has {} x {} values", getDataName(), _writeDataBuffer.gradient.rows(), _writeDataBuffer.gradient.cols());
  }
}

void WriteDataContext::storeBufferedData(double currentTime)
{
  _providedData->values()         = _writeDataBuffer.values;   // @todo this line should become unnecessary!
  _providedData->gradientValues() = _writeDataBuffer.gradient; // @todo this line should become unnecessary!
  _providedData->timeStepsStorage().setSampleAtTime(currentTime, _writeDataBuffer);
}

void WriteDataContext::clearStorage()
{
  // need to not only clear _providedData->timeStepsStorage(), but also _toData->timeStepsStorage() as soon as we map from storage to storage.
  _providedData->timeStepsStorage().clearAll();
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
