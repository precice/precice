#include "precice/impl/DataContext.hpp"
#include <memory>
#include "utils/EigenHelperFunctions.hpp"

namespace precice {
namespace impl {

logging::Logger DataContext::_log{"impl::DataContext"};

DataContext::DataContext(mesh::PtrData data, mesh::PtrMesh mesh)
{
  PRECICE_ASSERT(data);
  _providedData = data;
  PRECICE_ASSERT(mesh);
  _mesh = mesh;
}

std::string DataContext::getDataName() const
{
  PRECICE_ASSERT(_providedData);
  return _providedData->getName();
}

DataID DataContext::getFromDataID(int dataVectorIndex) const
{
  PRECICE_ASSERT(hasMapping());
  PRECICE_ASSERT(dataVectorIndex < _fromData.size())
  PRECICE_ASSERT(_fromData[dataVectorIndex]);
  return _fromData[dataVectorIndex]->getID();
}

void DataContext::resetData()
{
  // See also https://github.com/precice/precice/issues/1156.
  _providedData->toZero();
  if (hasMapping()) {
    PRECICE_ASSERT(hasWriteMapping());
    std::for_each(_toData.begin(), _toData.end(), [](auto &data) { data->toZero(); });
  }
}

DataID DataContext::getToDataID(int dataVectorIndex) const
{
  PRECICE_ASSERT(hasMapping());
  PRECICE_ASSERT(dataVectorIndex < _toData.size())
  PRECICE_ASSERT(_toData[dataVectorIndex]);
  return _toData[dataVectorIndex]->getID();
}

DataID DataContext::getDataDimensions() const
{
  PRECICE_ASSERT(_providedData);
  return _providedData->getDimensions();
}

std::string DataContext::getMeshName() const
{
  PRECICE_ASSERT(_mesh);
  return _mesh->getName();
}

MeshID DataContext::getMeshID() const
{
  PRECICE_ASSERT(_mesh);
  return _mesh->getID();
}

void DataContext::addMapping(MappingContext mappingContext, mesh::PtrData fromData, mesh::PtrData toData)
{
  PRECICE_ASSERT(fromData);
  PRECICE_ASSERT(toData);
  _mappingContext.emplace_back(mappingContext);
  PRECICE_ASSERT(fromData == _providedData || toData == _providedData, "Either fromData or toData has to equal _providedData.");
  PRECICE_ASSERT(fromData->getName() == getDataName());
  _fromData.emplace_back(fromData);
  PRECICE_ASSERT(toData->getName() == getDataName());
  _toData.emplace_back(toData);
  PRECICE_ASSERT(_toData != _fromData);
}

bool DataContext::hasMapping() const
{
  return hasReadMapping() || hasWriteMapping();
}

bool DataContext::isMappingRequired()
{
  if (not hasMapping()) {
    return false;
  }

  return std::any_of(_mappingContext.begin(), _mappingContext.end(), [](const auto &context) {
                const auto timing    = context.timing;
                const bool mapNow    = (timing == mapping::MappingConfiguration::ON_ADVANCE) || (timing == mapping::MappingConfiguration::INITIAL);
                return (mapNow && !context.hasMappedData); });
}

void DataContext::mapData()
{
  PRECICE_ASSERT(hasMapping());
  // Reset the toData
  std::for_each(_toData.begin(), _toData.end(), [](auto &data) { data->toZero(); });

  // Execute the mapping
  for (unsigned int i = 0; i < _mappingContext.size(); ++i) {
    const DataID fromDataID = getFromDataID(i);
    const DataID toDataID   = getToDataID(i);
    _mappingContext[i].mapping->map(fromDataID, toDataID);
    PRECICE_DEBUG("Mapped values = {}", utils::previewRange(3, _toData[i]->values()));
  }
}

bool DataContext::hasReadMapping() const
{
  return std::any_of(_toData.begin(), _toData.end(), [this](auto &data) { return data == _providedData; });
}

bool DataContext::hasWriteMapping() const
{
  return std::any_of(_fromData.begin(), _fromData.end(), [this](auto &data) { return data == _providedData; });
}

} // namespace impl
} // namespace precice
