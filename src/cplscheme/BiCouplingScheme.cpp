#include "BiCouplingScheme.hpp"
#include <algorithm>
#include <map>
#include <memory>
#include <ostream>
#include <type_traits>
#include <utility>
#include "cplscheme/BaseCouplingScheme.hpp"
#include "cplscheme/CouplingData.hpp"
#include "cplscheme/SharedPointer.hpp"
#include "logging/LogMacros.hpp"
#include "m2n/M2N.hpp"
#include "m2n/SharedPointer.hpp"
#include "mesh/Data.hpp"
#include "utils/Helpers.hpp"

namespace precice {
namespace cplscheme {

BiCouplingScheme::BiCouplingScheme(
    double                        maxTime,
    int                           maxTimeWindows,
    double                        timeWindowSize,
    int                           validDigits,
    const std::string &           firstParticipant,
    const std::string &           secondParticipant,
    const std::string &           localParticipant,
    m2n::PtrM2N                   m2n,
    int                           maxIterations,
    CouplingMode                  cplMode,
    constants::TimesteppingMethod dtMethod)
    : BaseCouplingScheme(maxTime, maxTimeWindows, timeWindowSize, validDigits, localParticipant, maxIterations, cplMode, dtMethod),
      _m2n(m2n),
      _firstParticipant(firstParticipant),
      _secondParticipant(secondParticipant)
{
  PRECICE_ASSERT(_firstParticipant != _secondParticipant,
                 "First participant and second participant must have different names.");
  if (localParticipant == _firstParticipant) {
    setDoesFirstStep(true);
  } else if (localParticipant == _secondParticipant) {
    setDoesFirstStep(false);
  } else {
    PRECICE_ERROR("Name of local participant \""
                  << localParticipant << "\" does not match any "
                  << "participant specified for the coupling scheme.");
  }
}

void BiCouplingScheme::addDataToSend(
    mesh::PtrData data,
    mesh::PtrMesh mesh,
    bool          requiresInitialization)
{
  PRECICE_TRACE();
  int id = data->getID();
  if (!utils::contained(id, _sendData)) {
    PtrCouplingData     ptrCplData(new CouplingData(data, mesh, requiresInitialization));
    DataMap::value_type pair = std::make_pair(id, ptrCplData);
    _sendData.insert(pair);
  } else {
    PRECICE_ERROR("Data \"" << data->getName() << "\" cannot be added twice for sending. Please remove any duplicate <exchange data=\"" << data->getName() << "\" .../> tags");
  }
}

void BiCouplingScheme::addDataToReceive(
    mesh::PtrData data,
    mesh::PtrMesh mesh,
    bool          requiresInitialization)
{
  PRECICE_TRACE();
  int id = data->getID();
  if (!utils::contained(id, _receiveData)) {
    PtrCouplingData     ptrCplData(new CouplingData(data, mesh, requiresInitialization));
    DataMap::value_type pair = std::make_pair(id, ptrCplData);
    _receiveData.insert(pair);
  } else {
    PRECICE_ERROR("Data \"" << data->getName() << "\" cannot be added twice for receiving. Please remove any duplicate <exchange data=\"" << data->getName() << "\" ... /> tags");
  }
}

std::vector<std::string> BiCouplingScheme::getCouplingPartners() const
{
  std::vector<std::string> partnerNames;
  // Add non-local participant
  if (doesFirstStep()) {
    partnerNames.push_back(_secondParticipant);
  } else {
    partnerNames.push_back(_firstParticipant);
  }
  return partnerNames;
}

bool BiCouplingScheme::receiveConvergence()
{
  PRECICE_ASSERT(doesFirstStep(), "For convergence information the receiving participant is always the first one.");
  bool convergence;
  _m2n->receive(convergence);
  return convergence;
}

CouplingData *BiCouplingScheme::getSendData(
    int dataID)
{
  PRECICE_TRACE(dataID);
  DataMap::iterator iter = _sendData.find(dataID);
  if (iter != _sendData.end()) {
    return &(*(iter->second));
  }
  return nullptr;
}

CouplingData *BiCouplingScheme::getReceiveData(
    int dataID)
{
  PRECICE_TRACE(dataID);
  DataMap::iterator iter = _receiveData.find(dataID);
  if (iter != _receiveData.end()) {
    return &(*(iter->second));
  }
  return nullptr;
}

} // namespace cplscheme
} // namespace precice