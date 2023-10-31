#include "WriteDataContext.hpp"

#include "utils/EigenHelperFunctions.hpp"

namespace precice::impl {

logging::Logger WriteDataContext::_log{"impl::WriteDataContext"};

WriteDataContext::WriteDataContext(mesh::PtrData data,
                                   mesh::PtrMesh mesh)
    : DataContext(data, mesh),
      _writeDataBuffer(data->getDimensions())
{
}

void WriteDataContext::resetData(bool atEndOfWindow, bool isTimeWindowComplete)
{
  PRECICE_ASSERT(atEndOfWindow || !isTimeWindowComplete, "isTimeWindowComplete without atEndOfWindow doesn't make sense!");

  // See also https://github.com/precice/precice/issues/1156.
  _providedData->toZero();

  // reset writeDataBuffer
  _writeDataBuffer.values.setZero();
  _writeDataBuffer.gradients.setZero();

  PRECICE_ASSERT(!hasReadMapping(), "Read mapping is not allowed for WriteDataContext.");
  // reset all toData
  if (hasWriteMapping()) {
    std::for_each(_mappingContexts.begin(), _mappingContexts.end(), [](auto &context) { context.toData->toZero(); });
  }

  if (!atEndOfWindow) {
    return;
  }

  if (isTimeWindowComplete && hasWriteMapping()) { // manually overwrite value at beginning with value from end. Need this exception for WriteDataContext with write mapping, because CouplingScheme is not able to update _providedData, if write mapping sits between _providedData and _toData. CouplingScheme in this case only has access to _toData.
    _providedData->timeStepsStorage().move();
    _providedData->timeStepsStorage().trim();
  } else {
    _providedData->timeStepsStorage().trim();
  }
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
    const auto     actualColumnSize   = static_cast<SizeType>(_writeDataBuffer.gradients.cols());

    // Shrink Buffer
    if (expectedColumnSize < actualColumnSize) {
      _writeDataBuffer.gradients.resize(spaceDimensions, expectedColumnSize);
    }

    // Enlarge Buffer
    if (expectedColumnSize > actualColumnSize) {
      const auto columnLeftToAllocate = expectedColumnSize - actualColumnSize;
      utils::append(_writeDataBuffer.gradients, Eigen::MatrixXd(Eigen::MatrixXd::Zero(spaceDimensions, columnLeftToAllocate)));
    }
    PRECICE_DEBUG("Gradient Data {} now has {} x {} values", getDataName(), _writeDataBuffer.gradients.rows(), _writeDataBuffer.gradients.cols());
  }
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
