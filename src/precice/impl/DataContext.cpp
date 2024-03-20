#include <iterator>
#include <memory>
#include <utility>

#include "precice/impl/DataContext.hpp"
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
  for (auto &kv : _initialGuesses) {
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
  return _mesh->nVertices();
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

int DataContext::mapData(std::optional<double> after, bool skipZero)
{
  PRECICE_TRACE(getMeshName(), getDataName());
  PRECICE_ASSERT(hasMapping());

  int executedMappings{0};

  // Execute the mappings
  for (auto &context : _mappingContexts) {
    PRECICE_ASSERT(!context.fromData->stamples().empty(),
                   "There must be samples at this point!");

    // linear lookup should be sufficient here
    const auto timestampExists = [times = context.toData->timeStepsStorage().getTimes()](double lookup) -> bool {
      return std::any_of(times.data(), std::next(times.data(), times.size()), [lookup](double time) {
        return math::equals(time, lookup);
      });
    };

    auto &mapping = *context.mapping;

    const auto dataDims = context.fromData->getDimensions();

    for (const auto &stample : context.fromData->stamples()) {
      // skip stamples before given time
      if (after && math::smallerEquals(stample.timestamp, *after)) {
        PRECICE_DEBUG("Skipping stample t={} (not after {})", stample.timestamp, *after);
        continue;
      }
      // skip existing stamples
      if (timestampExists(stample.timestamp)) {
        PRECICE_DEBUG("Skipping stample t={} (exists)", stample.timestamp);
        continue;
      }

      time::Sample outSample{
          dataDims,
          Eigen::VectorXd::Zero(dataDims * mapping.getOutputMesh()->nVertices())};

      bool skipMapping = skipZero && stample.sample.values.isZero();

      PRECICE_INFO("Mapping \"{}\" for t={} from \"{}\" to \"{}\"{}",
                   getDataName(), stample.timestamp, mapping.getInputMesh()->getName(), mapping.getOutputMesh()->getName(),
                   (skipMapping ? " (skipped zero sample)" : ""));
      if (!skipMapping) {
        if (mapping.requiresInitialGuess()) {
          const FromToDataIDs key{context.fromData->getID(), context.toData->getID()};
          mapping.map(stample.sample, outSample.values, _initialGuesses[key]);
        } else {
          mapping.map(stample.sample, outSample.values);
        }
        PRECICE_DEBUG("Mapped values (t={}) = {}", stample.timestamp, utils::previewRange(3, outSample.values));
        ++executedMappings;
      }

      // Store data from mapping buffer in storage
      context.toData->setSampleAtTime(stample.timestamp, std::move(outSample));
    }
  }
  return executedMappings;
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
