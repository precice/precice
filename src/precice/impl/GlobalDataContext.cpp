#include "precice/impl/GlobalDataContext.hpp"
#include <memory>
#include "mesh/GlobalData.hpp" //TODO: Should compile without this since GlobalData is already forward declared in mesh/SharedPointer.hpp, but doesn't.
#include "time/Waveform.hpp"
#include "utils/EigenHelperFunctions.hpp"

namespace precice::impl {

logging::Logger GlobalDataContext::_log{"impl::GlobalDataContext"};

GlobalDataContext::GlobalDataContext(
    mesh::PtrGlobalData data)
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

mesh::PtrGlobalData GlobalDataContext::providedData() const
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

} // namespace precice::impl
