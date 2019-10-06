#include <limits>

#include "cplscheme/CouplingData.hpp"
#include "math/math.hpp"
#include "acceleration/HierarchicalAitkenAcceleration.hpp"
#include "utils/EigenHelperFunctions.hpp"
#include "utils/Helpers.hpp"

namespace precice
{
namespace acceleration
{

HierarchicalAitkenAcceleration::HierarchicalAitkenAcceleration(
    double           initialRelaxation,
    std::vector<int> dataIDs)
    : _initialRelaxation(initialRelaxation),
      _dataIDs(dataIDs)
{
  PRECICE_CHECK((_initialRelaxation > 0.0) && (_initialRelaxation <= 1.0),
        "Initial relaxation factor for aitken acceleration has to "
        << "be larger than zero and smaller or equal than one!");
}

void HierarchicalAitkenAcceleration::initialize(DataMap &cplData)
{
  PRECICE_TRACE();
  PRECICE_CHECK(utils::contained(*_dataIDs.begin(), cplData),
        "Data with ID " << *_dataIDs.begin() << " is not contained in data given at initialization!");
  size_t entries = cplData[*_dataIDs.begin()]->values->size(); // Add zero boundaries
  PRECICE_ASSERT((entries - 1) % 2 == 0);                           // entries has to be an odd number
  double          initializer = std::numeric_limits<double>::max();
  Eigen::VectorXd toAppend    = Eigen::VectorXd::Constant(entries, initializer);
  utils::append(_residual, toAppend);

  size_t entriesCurrentLevel = 1;
  size_t totalEntries        = 2;               // Boundary entries
  _aitkenFactors.push_back(_initialRelaxation); // Boundary entries
  while (totalEntries < entries) {
    _aitkenFactors.push_back(_initialRelaxation);
    totalEntries += entriesCurrentLevel;
    entriesCurrentLevel *= 2;
  }
  PRECICE_ASSERT(totalEntries == entries);
  //  INFO ( "HierarchicalAitkenAcceleration: level count = " << _aitkenFactors.size() );

  // Append column for old values if not done by coupling scheme yet
  for (DataMap::value_type &pair : cplData) {
    int cols = pair.second->oldValues.cols();
    if (cols < 1) {
      utils::append(pair.second->oldValues, (Eigen::VectorXd) Eigen::VectorXd::Zero(pair.second->values->size()));
    }
  }
}

void HierarchicalAitkenAcceleration::performAcceleration(
    DataMap &cplData)
{
  PRECICE_TRACE();
  typedef Eigen::VectorXd DataValues;

  // Compute aitken relaxation factor
  PRECICE_ASSERT(utils::contained(*_dataIDs.begin(), cplData));
  auto &values = *cplData[*_dataIDs.begin()]->values;

  // Attention: no passing by ref any more --> needs to be written back at the end of the function!
  Eigen::VectorXd oldValues = cplData[*_dataIDs.begin()]->oldValues.col(0);

  // Compute current residuals
  DataValues residual = values;
  residual -= oldValues;

  // Compute residual deltas and temporarily store it in _residuals
  DataValues residualDelta = _residual;
  residualDelta *= -1.0;
  residualDelta += residual;

  std::vector<double> nominators(_aitkenFactors.size(), 0.0);
  std::vector<double> denominators(_aitkenFactors.size(), 0.0);

  //  PRECICE_INFO( "" );
  //  PRECICE_INFO( "values = " << values );
  //  PRECICE_INFO( "oldValues = " << oldValues );
  //  PRECICE_INFO( "residual = " << residual );
  //  PRECICE_INFO( "_residual = " << _residual );
  //  PRECICE_INFO( "residualDelta = " << residualDelta );

  // Hierarchize entries
  size_t entries             = residual.size();
  size_t treatedEntries      = 2;
  size_t entriesCurrentLevel = std::pow(2.0, (int) (_aitkenFactors.size() - 2));
  for (size_t level = _aitkenFactors.size() - 1; level > 0; level--) {
    size_t stepsize = (entries - 1) / std::pow(2.0, (int) (level - 1));
    PRECICE_ASSERT(stepsize % 2 == 0);
    size_t index = stepsize / 2;

    for (size_t i = 0; i < entriesCurrentLevel; i++) {
      _residual(index) -= (_residual(index - stepsize / 2) +
                           _residual(index + stepsize / 2)) /
                          2.0;
      residualDelta(index) -= (residualDelta(index - stepsize / 2) +
                               residualDelta(index + stepsize / 2)) /
                              2.0;
      values(index) -= (values(index - stepsize / 2) +
                        values(index + stepsize / 2)) /
                       2.0;
      oldValues(index) -= (oldValues(index - stepsize / 2) +
                           oldValues(index + stepsize / 2)) /
                          2.0;
      index += stepsize;
    }
    treatedEntries += entriesCurrentLevel;
    entriesCurrentLevel /= 2;
  }
  PRECICE_ASSERT(treatedEntries == entries);
  //  PRECICE_INFO( "hierarchized values = " << values );
  //  PRECICE_INFO( "hierarchized oldValues = " << oldValues );
  //  PRECICE_INFO( "hierarchized _residual = " << _residual );
  //  PRECICE_INFO( "hierarchized residualDelta = " << residualDelta );

  // Compute and perform relaxation with aitken factor
  nominators[0] = _residual(0) * residualDelta(0) +
                  _residual(entries - 1) * residualDelta(entries - 1);
  denominators[0] = residualDelta(0) * residualDelta(0) +
                    residualDelta(entries - 1) * residualDelta(entries - 1);
  computeAitkenFactor(0, nominators[0], denominators[0]);
  double omega         = _aitkenFactors[0];
  double oneMinusOmega = 1.0 - omega;
  for (DataMap::value_type &pair : cplData) {
    auto &      values    = *pair.second->values;
    const auto &oldValues = pair.second->oldValues.col(0);
    values(0)             = values(0) * omega + oldValues(0) * oneMinusOmega;
    values(entries - 1)   = values(entries - 1) * omega + oldValues(entries - 1) * oneMinusOmega;
  }
  treatedEntries      = 2;
  entriesCurrentLevel = 1;
  for (size_t level = 1; level < _aitkenFactors.size(); level++) {
    size_t stepsize = (entries - 1) / std::pow(2.0, (int) (level - 1));
    size_t index    = stepsize / 2;
    for (size_t i = 0; i < entriesCurrentLevel; i++) {
      nominators[level] += _residual(index) * residualDelta(index);
      denominators[level] += residualDelta(index) * residualDelta(index);
      index += stepsize;
    }
    computeAitkenFactor(level, nominators[level], denominators[level]);
    omega         = _aitkenFactors[level];
    oneMinusOmega = 1.0 - omega;
    //for ( DataMap::value_type & pair : cplData ) {
    //  DataValues & values = *pair.second.values;
    //  DataValues & oldValues = pair.second.oldValues.getColumn(0);
    index = stepsize / 2;
    for (size_t i = 0; i < entriesCurrentLevel; i++) {
      values(index) = values(index) * omega + oldValues(index) * oneMinusOmega;
      index += stepsize;
    }
    //}
    treatedEntries += entriesCurrentLevel;
    entriesCurrentLevel *= 2;
  }
  PRECICE_ASSERT(treatedEntries == entries);

  _residual = residual; // Overwrite old residual by current one

  //  double nom = 0.0;
  //  double denom = 0.0;
  //  for ( size_t i=0; i < entries; i++ ) {
  //    nom += _residual[i] * residualDelta[i];
  //    denom += residualDelta[i] * residualDelta[i];
  //  }
  //  computeAitkenFactor ( 0, nom, denom );
  //  double omega = _aitkenFactors[0];
  //  double oneMinusOmega = 1.0 - omega;
  //  for ( DataMap::value_type & pair : cplData ) {
  //    DataValues & values = *pair.second.values;
  //    DataValues & oldValues = pair.second.oldValues.getColumn(0);
  //    for ( size_t i=0; i < entries; i++ ) {
  //      values[i] = values[i] * omega + oldValues[i] * oneMinusOmega;
  //    }
  //  }
  //  _residual = residual; // Overwrite old residual by current one

  //  PRECICE_INFO( "relaxed hierarchized values = " << values );

  // Dehierarchization
  treatedEntries      = 2;
  entriesCurrentLevel = 1;
  for (size_t level = 1; level < _aitkenFactors.size(); level++) {
    size_t stepsize = (entries - 1) / std::pow(2.0, (int) (level - 1));
    size_t index    = stepsize / 2;
    for (size_t i = 0; i < entriesCurrentLevel; i++) {
      values(index) +=
          (values(index - stepsize / 2) + values(index + stepsize / 2)) / 2.0;
      oldValues(index) +=
          (oldValues(index - stepsize / 2) + oldValues(index + stepsize / 2)) / 2.0;
      index += stepsize;
    }
    treatedEntries += entriesCurrentLevel;
    entriesCurrentLevel *= 2;
  }
  PRECICE_ASSERT(treatedEntries == entries);
  //  PRECICE_INFO( "relaxed values = " << values );
  //  PRECICE_INFO( "oldValues = " << oldValues );

  // save back oldValues in cplData. Eigen does not allow call by ref for blocks (cols, rows), thus explicitly write back
  cplData[*_dataIDs.begin()]->oldValues.col(0) = oldValues;

  _iterationCounter++;
}

void HierarchicalAitkenAcceleration::iterationsConverged(
    DataMap &cplData)
{
  _iterationCounter = 0;
  _residual         = Eigen::VectorXd::Constant(_residual.size(), std::numeric_limits<double>::max());
}

void HierarchicalAitkenAcceleration::computeAitkenFactor(
    size_t level,
    double nominator,
    double denominator)
{
  // Select/compute aitken factor depending on current iteration count
  if (_iterationCounter == 0) {
    //PRECICE_INFO( "First iteration (nom = " << nominator << ", den = " << denominator );
    _aitkenFactors[level] = math::sign(_aitkenFactors[level]) *
                            std::min(_initialRelaxation, std::abs(_aitkenFactors[level]));
  } else {
    //PRECICE_INFO( "Aitken factor = - " << _aitkenFactors[level] << " * "
    //               << nominator << " / " << denominator );
    if (math::equals(std::sqrt(denominator), 0.0)) {
      _aitkenFactors[level] = 1.0;
    } else {
      _aitkenFactors[level] = -_aitkenFactors[level] * (nominator / denominator);
    }
  }
  //PRECICE_INFO( "Level " << level << " relaxation factor = " << _aitkenFactors[level] );
}

/** ---------------------------------------------------------------------------------------------
 *         getDesignSpecification()
 *
 * @brief: Returns the design specification corresponding to the given coupling data.
 *         This information is needed for convergence measurements in the coupling scheme.
 *  ---------------------------------------------------------------------------------------------
 */
std::map<int, Eigen::VectorXd> HierarchicalAitkenAcceleration::getDesignSpecification(
    DataMap &cplData)
{
  PRECICE_ERROR("Design specification for Aitken relaxation is not supported yet.");

  std::map<int, Eigen::VectorXd> designSpecifications;
  int                            off = 0;
  for (int id : _dataIDs) {
    int             size = cplData[id]->values->size();
    Eigen::VectorXd q    = Eigen::VectorXd::Zero(size);
    for (int i = 0; i < size; i++) {
      q(i) = _designSpecification(i + off);
    }
    off += size;
    std::map<int, Eigen::VectorXd>::value_type pair = std::make_pair(id, q);
    designSpecifications.insert(pair);
  }
  return designSpecifications;
}

void HierarchicalAitkenAcceleration::setDesignSpecification(
    Eigen::VectorXd &q)
{
  _designSpecification = q;
  PRECICE_ERROR("design specification for Aitken relaxation is not supported yet.");
}
}
} // namespace precice, acceleration
