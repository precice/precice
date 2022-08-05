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

DataID DataContext::getFromDataID(size_t dataVectorIndex) const
{
  PRECICE_ASSERT(hasMapping());
  PRECICE_ASSERT(dataVectorIndex < _fromData.size())
  PRECICE_ASSERT(_fromData[dataVectorIndex]);
  return _fromData[dataVectorIndex]->getID();
}

DataID DataContext::getToDataID(size_t dataVectorIndex) const
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

void DataContext::appendMapping(MappingContext mappingContext, mesh::PtrData fromData, mesh::PtrData toData)
{
  PRECICE_ASSERT(fromData);
  PRECICE_ASSERT(toData);
  // Make sure we don't append a mapping twice
#ifndef NDEBUG
  for (unsigned int i = 0; i < _mappingContexts.size(); ++i) {
    PRECICE_ASSERT(!((_mappingContexts[i].mapping == mappingContext.mapping) && (_fromData[i] == fromData) && (_toData[i] == toData)), "The appended mapping already exists.");
  }
#endif
  _mappingContexts.emplace_back(mappingContext);
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

  PRECICE_ASSERT(std::all_of(_mappingContexts.begin(), _mappingContexts.end(), [this](const auto &context) { return context.timing == _mappingContexts[0].timing; }), "Different mapping timings for the same data context are not supported");

  return std::any_of(_mappingContexts.begin(), _mappingContexts.end(), [](const auto &context) {
    const auto timing = context.timing;
    const bool mapNow = (timing == mapping::MappingConfiguration::ON_ADVANCE) || (timing == mapping::MappingConfiguration::INITIAL);
    return (mapNow && !context.hasMappedData); });
}

void DataContext::mapData()
{
  PRECICE_ASSERT(hasMapping());
  // Execute the mapping
  for (unsigned int i = 0; i < _mappingContexts.size(); ++i) {
    const DataID fromDataID = getFromDataID(i);
    const DataID toDataID   = getToDataID(i);
    // Reset the toData before executing the mapping
    _toData[i]->toZero();
    _mappingContexts[i].mapping->map(fromDataID, toDataID);
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
