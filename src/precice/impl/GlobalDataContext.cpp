#include "precice/impl/GlobalDataContext.hpp"
#include <memory>
// #include "mesh/GlobalData.hpp" //TODO: Should compile without this since GlobalData is already forward declared in mesh/SharedPointer.hpp, but doesn't.
#include "mesh/Data.hpp"
#include "time/Waveform.hpp"
#include "utils/EigenHelperFunctions.hpp"

namespace precice::impl {

logging::Logger GlobalDataContext::_log{"impl::GlobalDataContext"};

GlobalDataContext::GlobalDataContext(
    mesh::PtrData data,
    std::string   direction)
{
  PRECICE_ASSERT(data);
  _providedData = data;
}

std::string GlobalDataContext::getDataName() const
{
  PRECICE_ASSERT(_providedData);
  return _providedData->getName();
}

DataID GlobalDataContext::getDataDimensions() const
{
  PRECICE_ASSERT(_providedData);
  return _providedData->getDimensions();
}

void GlobalDataContext::resetData()
{
  // See also https://github.com/precice/precice/issues/1156.
  _providedData->toZero();
}

mesh::PtrData GlobalDataContext::providedData() const
{
  PRECICE_ASSERT(_providedData);
  return _providedData;
}

int GlobalDataContext::getInterpolationOrder() const
{
  return _waveform->getInterpolationOrder();
}

Eigen::VectorXd GlobalDataContext::sampleWaveformAt(double normalizedDt)
{
  return _waveform->sample(normalizedDt);
}

void GlobalDataContext::initializeWaveform()
{
  _waveform->initialize(_providedData->values());
}

void GlobalDataContext::moveToNextWindow()
{
  _waveform->moveToNextWindow();
}

std::string GlobalDataContext::getDirection()
{
  return _direction;
}

void GlobalDataContext::storeDataInWaveform()
{
  _waveform->store(_providedData->values()); // store mapped or received _providedData in the _waveform
}

} // namespace precice::impl
