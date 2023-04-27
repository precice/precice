#include "ReadDataContext.hpp"

#include "time/Waveform.hpp"

namespace precice::impl {

logging::Logger ReadDataContext::_log{"impl::ReadDataContext"};

ReadDataContext::ReadDataContext(
    mesh::PtrData data,
    mesh::PtrMesh mesh,
    int           interpolationOrder)
    : DataContext(data, mesh)
{
  _waveform = std::make_shared<time::Waveform>(interpolationOrder);
}

void ReadDataContext::appendMappingConfiguration(MappingContext &mappingContext, const MeshContext &meshContext)
{
  PRECICE_ASSERT(!hasReadMapping(), "The read data context must be unique. Otherwise we would have an ambiguous read data operation on the user side.")
  PRECICE_ASSERT(meshContext.mesh->hasDataName(getDataName()));
  mesh::PtrData data = meshContext.mesh->data(getDataName());
  PRECICE_ASSERT(data != _providedData, "Data the read mapping is mapping from needs to be different from _providedData");
  mappingContext.fromData = data;
  mappingContext.toData   = _providedData;
  appendMapping(mappingContext);
  PRECICE_ASSERT(hasReadMapping());
}

Eigen::VectorXd ReadDataContext::readValues(std::vector<int> indices, double normalizedDt)
{
  const auto vertexCount    = getDataSize() / getDataDimensions();
  const auto valuesInternal = _waveform->sample(normalizedDt);
  auto       values         = Eigen::VectorXd(indices.size() * getDataDimensions());
  for (int i = 0; i < indices.size(); i++) {
    const auto valueIndex = indices[i];
    PRECICE_CHECK(0 <= valueIndex && valueIndex < vertexCount,
                  "Cannot read data \"{}\" from invalid Vertex ID ({}). "
                  "Please make sure you only use the results from calls to setMeshVertex/Vertices().",
                  getDataName(), valueIndex);
    int offsetInternal = valueIndex * getDataDimensions();
    int offset         = i * getDataDimensions();
    for (int dim = 0; dim < getDataDimensions(); dim++) {
      values[offset + dim] = valuesInternal[offsetInternal + dim];
    }
  }
  return values;
}

int ReadDataContext::getInterpolationOrder() const
{
  return _waveform->getInterpolationOrder();
}

void ReadDataContext::storeDataInWaveform()
{
  _waveform->store(_providedData->values()); // store mapped or received _providedData in the _waveform
}

void ReadDataContext::initializeWaveform()
{
  PRECICE_ASSERT(not hasWriteMapping(), "Write mapping does not need waveforms.");
  _waveform->initialize(_providedData->values());
}

void ReadDataContext::moveToNextWindow()
{
  _waveform->moveToNextWindow();
}

} // namespace precice::impl
