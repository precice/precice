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

void WriteDataContext::writeValuesIntoDataBuffer(::precice::span<const VertexID> vertices, ::precice::span<const double> values)
{
  Eigen::Map<const Eigen::MatrixXd> inputData(values.data(), getDataDimensions(), vertices.size());
  Eigen::Map<Eigen::MatrixXd>       localData(_writeDataBuffer.values.data(), getDataDimensions(), getMesh().vertices().size());

  for (int i = 0; i < vertices.size(); ++i) {
    const auto vid = vertices[i];
    PRECICE_CHECK(getMesh().isValidVertexID(vid),
                  "Cannot write data \"{}\" to invalid Vertex ID ({}) of mesh \"{}\". Please make sure you only use the results from calls to setMeshVertex/Vertices().",
                  getDataName(), vid, getMeshName());
    localData.col(vid) = inputData.col(i);
  }
}

void WriteDataContext::writeGradientsIntoDataBuffer(::precice::span<const VertexID> vertices, ::precice::span<const double> gradients)
{
  const auto                        gradientComponents = getMesh().getDimensions() * getDataDimensions();
  Eigen::Map<const Eigen::MatrixXd> inputGradients(gradients.data(), gradientComponents, vertices.size());
  Eigen::Map<Eigen::MatrixXd>       localGradients(_writeDataBuffer.gradients.data(), gradientComponents, getMesh().vertices().size());

  for (int i = 0; i < vertices.size(); ++i) {
    const auto vid = vertices[i];
    PRECICE_CHECK(getMesh().isValidVertexID(vid),
                  "Cannot write gradient for data \"{}\" to invalid Vertex ID ({}) of mesh \"{}\". Please make sure you only use the results from calls to setMeshVertex/Vertices().",
                  getDataName(), vid, getMeshName());
    localGradients.col(vid) = inputGradients.col(i);
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
    const SizeType spaceDimensions = getMesh().getDimensions();

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
  PRECICE_DEBUG("storeBufferedData {} at time {}", getDataName(), currentTime);
  _providedData->setSampleAtTime(currentTime, _writeDataBuffer);
}

void WriteDataContext::clearStorage()
{
  // need to not only clear _providedData->timeStepsStorage(), but also _toData->timeStepsStorage() as soon as we map from storage to storage.
  _providedData->timeStepsStorage().clearAll();
}

mesh::PtrData WriteDataContext::providedData()
{
  PRECICE_ASSERT(_providedData);
  return _providedData;
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
