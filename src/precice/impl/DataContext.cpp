#include "precice/impl/DataContext.hpp"
#include <memory>
#include <utility>
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

void DataContext::resetInitialGuesses()
{
  for (auto &kv : _lastSolutions) {
    kv.second.setZero();
  }
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

std::string DataContext::getMeshName() const
{
  PRECICE_ASSERT(_mesh);
  return _mesh->getName();
}

int DataContext::getMeshVertexCount() const
{
  PRECICE_ASSERT(_mesh);
  return _mesh->vertices().size();
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
    context.clearToDataStorage();

    PRECICE_ASSERT(context.fromData->stamples().size() > 0);
    for (auto &stample : context.fromData->stamples()) {
      // Put data from storage into mapping buffer
      context.fromData->sample() = stample.sample;

      // Reset the toData before executing the mapping
      context.toData->toZero();
      const DataID fromDataID = context.fromData->getID();
      const DataID toDataID   = context.toData->getID();
      auto &       mapping    = *context.mapping;

      if (mapping.isTransient()) {
        auto key = std::make_pair(fromDataID, toDataID);
        if (_lastSolutions.count(key) == 0) {
          _lastSolutions.emplace(key, Eigen::VectorXd{});
        }
        auto &lastSolution = _lastSolutions[key];
        mapping.map(fromDataID, toDataID, lastSolution);
      } else {
        mapping.map(fromDataID, toDataID);
      }

      // Store data from mapping buffer in storage
      context.toData->setSampleAtTime(stample.timestamp, context.toData->sample());

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

bool DataContext::isValidVertexID(const VertexID id) const
{
  PRECICE_ASSERT(_mesh);
  return _mesh->isValidVertexID(id);
}

} // namespace precice::impl
