#include "Data.hpp"
#include <algorithm>
#include <utility>

#include "math/differences.hpp"
#include "precice/impl/Types.hpp"
#include "utils/EigenHelperFunctions.hpp"
#include "utils/assertion.hpp"

namespace precice::mesh {

Data::Data(
    std::string                        name,
    DataID                             id,
    int                                dimensions,
    int                                spatialDimensions,
    int                                waveformDegree,
    std::vector<std::optional<double>> lowerBound,
    std::vector<std::optional<double>> upperBound)
    : _waveform(waveformDegree),
      _lowerBound(lowerBound),
      _upperBound(upperBound),
      _name(std::move(name)),
      _id(id),
      _dimensions(dimensions),
      _spatialDimensions(spatialDimensions),
      _sample(_dimensions)
{
  PRECICE_ASSERT(dimensions > 0, dimensions);
}

Eigen::VectorXd &Data::values()
{
  return _sample.values;
}

const Eigen::VectorXd &Data::values() const
{
  return _sample.values;
}

const Eigen::MatrixXd &Data::gradients() const
{
  return _sample.gradients;
}

const time::Sample &Data::sample() const
{
  return _sample;
}

time::SampleResult Data::sampleAtTime(double time) const
{
  return _waveform.sample(time);
}

int Data::getWaveformDegree() const
{
  return _waveform.getInterpolationDegree();
}

time::Waveform &Data::waveform()
{
  return _waveform;
}

std::vector<std::optional<double>> Data::getLowerBound() const
{
  return _lowerBound;
}

std::vector<std::optional<double>> Data::getUpperBound() const
{
  return _upperBound;
}

void Data::moveToNextWindow()
{
  if (stamples().size() > 1) { // Needed to avoid CompositionalCouplingScheme callong moveToNextWindow on same Data multiple times. Could be simplified by replacing Waveform::move() with clearBefore(double time). See https://github.com/precice/precice/issues/1821.
    waveform().move();
    PRECICE_ASSERT(stamples().size() == 1);
    setGlobalSample(stamples().back().sample); // @todo at some point we should not need this anymore, when mapping, acceleration ... directly work on _timeStepsWaveform, see https://github.com/precice/precice/issues/1645
  }
}

void Data::setSampleAtTime(double time, const time::Sample &sample)
{
  PRECICE_ASSERT(sample.dataDims == getDimensions(), "Sample has incorrect data dimension");
  setGlobalSample(sample); // @todo at some point we should not need this anymore, when mapping, acceleration ... directly work on _timeStepsWaveform, see https://github.com/precice/precice/issues/1645
  _waveform.setSampleAtTime(time, sample);
}

void Data::setGlobalSample(const time::Sample &sample)
{
  PRECICE_ASSERT(not sample.values.hasNaN());
  _sample = sample;
}

void Data::emplaceSampleAtTime(double time)
{
  setSampleAtTime(time, time::Sample{getDimensions()});
}

void Data::emplaceSampleAtTime(double time, std::initializer_list<double> values)
{
  setSampleAtTime(time, time::Sample{getDimensions(),
                                     Eigen::Map<const Eigen::VectorXd>(values.begin(), values.size())});
}

void Data::emplaceSampleAtTime(double time, std::initializer_list<double> values, std::initializer_list<double> gradients)
{
  PRECICE_ASSERT(gradients.size() == values.size() * getSpatialDimensions(), "Gradient isn't correctly sized", values.size(), gradients.size());
  auto nVertices = values.size() / getDimensions();
  setSampleAtTime(time, time::Sample{getDimensions(),
                                     Eigen::Map<const Eigen::VectorXd>(values.begin(), values.size()),
                                     Eigen::Map<const Eigen::MatrixXd>(gradients.begin(), getSpatialDimensions(), nVertices * getDimensions())});
}

const std::string &Data::getName() const
{
  return _name;
}

DataID Data::getID() const
{
  return _id;
}

bool Data::hasGradient() const
{
  return _hasGradient;
}

bool Data::hasSamples() const
{
  return !_waveform.empty();
}

void Data::requireDataGradient()
{
  _hasGradient = true;
};

int Data::getDimensions() const
{
  return _dimensions;
}

void Data::allocateValues(int expectedCount)
{
  using SizeType = std::remove_cv<decltype(expectedCount)>::type;
  // Allocate data values
  const SizeType expectedSize = expectedCount * _dimensions;
  const auto     actualSize   = static_cast<SizeType>(_sample.values.size());
  // Shrink Buffer
  if (expectedSize < actualSize) {
    _sample.values.resize(expectedSize);
  }
  // Enlarge Buffer
  if (expectedSize > actualSize) {
    const auto leftToAllocate = expectedSize - actualSize;
    utils::append(_sample.values, Eigen::VectorXd(Eigen::VectorXd::Zero(leftToAllocate)));
  }
  PRECICE_DEBUG("Data {} now has {} values", _name, _sample.values.size());

  // Allocate gradient data values
  if (_hasGradient) {
    const SizeType spaceDimensions = _spatialDimensions;

    const SizeType expectedColumnSize = expectedCount * _dimensions;
    const auto     actualColumnSize   = static_cast<SizeType>(_sample.gradients.cols());

    // Shrink Buffer
    if (expectedColumnSize < actualColumnSize) {
      _sample.gradients.resize(spaceDimensions, expectedColumnSize);
    }

    // Enlarge Buffer
    if (expectedColumnSize > actualColumnSize) {
      const auto columnLeftToAllocate = expectedColumnSize - actualColumnSize;
      utils::append(_sample.gradients, Eigen::MatrixXd(Eigen::MatrixXd::Zero(spaceDimensions, columnLeftToAllocate)));
    }
    PRECICE_DEBUG("Gradient Data {} now has {} x {} values", _name, _sample.gradients.rows(), _sample.gradients.cols());
  }
}

int Data::getSpatialDimensions() const
{
  return _spatialDimensions;
}

} // namespace precice::mesh
