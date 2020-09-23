#include "acceleration/ConstantRelaxationAcceleration.hpp"
#include <Eigen/Core>
#include <map>
#include <memory>
#include <ostream>
#include "cplscheme/CouplingData.hpp"
#include "logging/LogMacros.hpp"
#include "utils/EigenHelperFunctions.hpp"
#include "utils/Helpers.hpp"
#include "utils/assertion.hpp"

namespace precice {
namespace acceleration {

ConstantRelaxationAcceleration::ConstantRelaxationAcceleration(
    double           relaxation,
    std::vector<int> dataIDs)
    : _relaxation(relaxation),
      _dataIDs(dataIDs)
{
  PRECICE_CHECK((relaxation > 0.0) && (relaxation <= 1.0),
                "Relaxation factor for constant relaxation acceleration "
                    << "has to be larger than zero and smaller or equal to one. Current relaxation factor is: " << relaxation);
}

void ConstantRelaxationAcceleration::initialize(DataMap &cplData)
{
  checkDataIDs(cplData);

  // Append column for old values if not done by coupling scheme yet
  int entries = 0;
  for (auto &elem : _dataIDs) {
    entries += cplData[elem]->values().size();
  }
  for (DataMap::value_type &pair : cplData) {
    int cols = pair.second->oldValues.cols();
    if (cols < 1) {
      PRECICE_ASSERT(pair.second->values().size() > 0, pair.first);
      utils::append(pair.second->oldValues, (Eigen::VectorXd) Eigen::VectorXd::Zero(pair.second->values().size()));
    }
  }
}

void ConstantRelaxationAcceleration::performAcceleration(DataMap &cplData)
{
  PRECICE_TRACE();
  double omega         = _relaxation;
  double oneMinusOmega = 1.0 - omega;
  for (DataMap::value_type &pair : cplData) {
    auto &      values    = pair.second->values();
    const auto &oldValues = pair.second->oldValues.col(0);
    values *= omega;
    values += oldValues * oneMinusOmega;
    PRECICE_DEBUG("pp values" << values);
  }
}

} // namespace acceleration
} // namespace precice
