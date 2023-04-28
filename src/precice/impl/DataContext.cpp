#include "precice/impl/DataContext.hpp"
#include <memory>
#include "utils/EigenHelperFunctions.hpp"

namespace precice::impl {

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

void DataContext::resetData()
{
  // See also https://github.com/precice/precice/issues/1156.
  _providedData->toZero();
  if (hasMapping()) {
    PRECICE_ASSERT(hasWriteMapping());
    PRECICE_ASSERT(!hasReadMapping());
    std::for_each(_mappingContexts.begin(), _mappingContexts.end(), [](auto &context) { context.toData->toZero(); });
  }
  // @todo need to also reset/clear (write)DataBuffer?
}

int DataContext::getDataDimensions() const
{
  PRECICE_ASSERT(_providedData);
  return _providedData->getDimensions();
}

int DataContext::getSpatialDimensions() const
{
  PRECICE_ASSERT(_providedData);
  return _providedData->getSpatialDimensions();
}

int DataContext::getDataSize() const
{
  PRECICE_ASSERT(_providedData);
  return _providedData->values().size();
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

bool DataContext::hasGradient() const
{
  PRECICE_ASSERT(_providedData);
  return _providedData->hasGradient();
}

void DataContext::appendMapping(MappingContext mappingContext)
{
  PRECICE_ASSERT(mappingContext.fromData);
  PRECICE_ASSERT(mappingContext.toData);
  // Make sure we don't append a mapping twice
#ifndef NDEBUG
  for (auto &context : _mappingContexts) {
    PRECICE_ASSERT(!((context.mapping == mappingContext.mapping) && (context.fromData == mappingContext.fromData) && (context.fromData == mappingContext.toData)), "The appended mapping already exists.");
  }
#endif
  _mappingContexts.emplace_back(mappingContext);
  PRECICE_ASSERT(mappingContext.fromData == _providedData || mappingContext.toData == _providedData, "Either fromData or toData has to equal _providedData.");
  PRECICE_ASSERT(mappingContext.fromData->getName() == getDataName());
  PRECICE_ASSERT(mappingContext.toData->getName() == getDataName());
}

bool DataContext::hasMapping() const
{
  return hasReadMapping() || hasWriteMapping();
}

void DataContext::mapData()
{
  PRECICE_ASSERT(hasMapping());
  // Execute the mapping
  for (auto &context : _mappingContexts) {
    // @todo messy. Try to improve this. Current problem: With clear all we also remove the data at WINDOW_START, which is not received by the coupling scheme.
    if (context.toData->timeStepsStorage().nTimes() > 0) {
      if (context.toData->timeStepsStorage().getTimes()[0] != time::Storage::WINDOW_START) {
        context.toData->timeStepsStorage().clearAll();
      } else {
        context.toData->timeStepsStorage().clear();
      }
    }

    PRECICE_ASSERT(context.fromData->getStamples().size() > 0);
    for (auto &stample : context.fromData->getStamples()) {
      // Put data from storage into mapping buffer
      context.fromData->sample() = stample.sample;

      // Reset the toData before executing the mapping
      context.toData->toZero();
      const DataID fromDataID = context.fromData->getID();
      const DataID toDataID   = context.toData->getID();
      context.mapping->map(fromDataID, toDataID);

      // Store data from mapping buffer in storage
      context.toData->timeStepsStorage().setSampleAtTime(stample.timestamp, context.toData->sample());

      PRECICE_DEBUG("Mapped values = {}", utils::previewRange(3, context.toData->values()));
    }
  }
}

bool DataContext::hasReadMapping() const
{
  return std::any_of(_mappingContexts.begin(), _mappingContexts.end(), [this](auto &context) { return context.toData == _providedData; });
}

bool DataContext::hasWriteMapping() const
{
  return std::any_of(_mappingContexts.begin(), _mappingContexts.end(), [this](auto &context) { return context.fromData == _providedData; });
}

} // namespace precice::impl
