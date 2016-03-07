// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "HierarchicalAitkenPostProcessing.hpp"
#include "../CouplingData.hpp"
#include "utils/Globals.hpp"
#include "utils/EigenHelperFunctions.hpp"
#include <limits>

namespace precice {
namespace cplscheme {
namespace impl {

tarch::logging::Log HierarchicalAitkenPostProcessing::
  _log ( "precice::cplscheme::HierarchicalAitkenPostProcessing" );

HierarchicalAitkenPostProcessing:: HierarchicalAitkenPostProcessing
(
  double initialRelaxation,
  std::vector<int>    dataIDs )
:
  _initialRelaxation ( initialRelaxation ),
  _dataIDs ( dataIDs ),
  _aitkenFactors (),
  _iterationCounter (),
  _residual (),
  _designSpecification()
{
  preciceCheck (
    (_initialRelaxation > 0.0) && (_initialRelaxation <= 1.0),
    "HierarchicalAitkenPostProcessing()",
    "Initial relaxation factor for aitken post processing has to "
    << "be larger than zero and smaller or equal than one!" );
}

void HierarchicalAitkenPostProcessing:: initialize
(
  DataMap & cplData )
{
  preciceTrace ( "initialize()" );
  preciceCheck ( utils::contained(*_dataIDs.begin(), cplData), "initialize()",
                 "Data with ID " << *_dataIDs.begin()
                 << " is not contained in data given at initialization!" );
  size_t entries = cplData[*_dataIDs.begin()]->values->size(); // Add zero boundaries
  assertion ( (entries - 1) % 2 == 0  ); // entries has to be an odd number
  double initializer = std::numeric_limits<double>::max ();
  Eigen::VectorXd toAppend = Eigen::VectorXd::Constant(entries, initializer);
  utils::append(_residual, toAppend);

  size_t entriesCurrentLevel = 1;
  size_t totalEntries = 2; // Boundary entries
  _aitkenFactors.push_back ( _initialRelaxation ); // Boundary entries
  while ( totalEntries < entries ) {
    _aitkenFactors.push_back ( _initialRelaxation );
    totalEntries += entriesCurrentLevel;
    entriesCurrentLevel *= 2;
  }
  assertion ( totalEntries == entries );
//  precicePrint ( "HierarchicalAitkenPostProcessing: level count = " << _aitkenFactors.size() );

  // Append column for old values if not done by coupling scheme yet
  for (DataMap::value_type& pair : cplData) {
    int cols = pair.second->oldValues.cols();
    if (cols < 1){
      utils::append(pair.second->oldValues, (Eigen::VectorXd) Eigen::VectorXd::Zero(pair.second->values->size()));
    }
  }
}

void HierarchicalAitkenPostProcessing:: performPostProcessing
(
  DataMap & cplData )
{
  preciceTrace ( "performPostProcessing()" );
  typedef Eigen::VectorXd DataValues;

  // Compute aitken relaxation factor
  assertion ( utils::contained(*_dataIDs.begin(), cplData) );
  auto& values = *cplData[*_dataIDs.begin()]->values;

  // Attention: no passing by ref any more --> needs to be written back at the end of the function!
  Eigen::VectorXd oldValues = cplData[*_dataIDs.begin()]->oldValues.col(0);

  // Compute current residuals
  DataValues residual = values;
  residual -= oldValues;

  // Compute residual deltas and temporarily store it in _residuals
  DataValues residualDelta = _residual;
  residualDelta *= -1.0;
  residualDelta += residual;

  std::vector<double> nominators ( _aitkenFactors.size(), 0.0 );
  std::vector<double> denominators ( _aitkenFactors.size(), 0.0 );

//  precicePrint ( "" );
//  precicePrint ( "values = " << values );
//  precicePrint ( "oldValues = " << oldValues );
//  precicePrint ( "residual = " << residual );
//  precicePrint ( "_residual = " << _residual );
//  precicePrint ( "residualDelta = " << residualDelta );

  // Hierarchize entries
  size_t entries = residual.size ();
  size_t treatedEntries = 2;
  size_t entriesCurrentLevel = std::pow(2.0, (int)(_aitkenFactors.size() - 2));
  for ( size_t level=_aitkenFactors.size()-1; level > 0; level-- ) {
    size_t stepsize = (entries - 1) / std::pow(2.0, (int)(level-1));
    assertion ( stepsize % 2 == 0 );
    size_t index = stepsize / 2;

    for ( size_t i=0; i < entriesCurrentLevel; i++ ) {
      _residual(index) -= ( _residual(index - stepsize/2) +
                            _residual(index + stepsize/2) ) / 2.0;
      residualDelta(index) -= ( residualDelta(index - stepsize/2) +
                                residualDelta(index + stepsize/2) ) / 2.0;
      values(index) -= ( values(index - stepsize/2) +
                         values(index + stepsize/2) ) / 2.0;
      oldValues(index) -= ( oldValues(index - stepsize/2) +
                            oldValues(index + stepsize/2) ) / 2.0;
      index += stepsize;
    }
    treatedEntries += entriesCurrentLevel;
    entriesCurrentLevel /= 2;
  }
  assertion ( treatedEntries == entries );
//  precicePrint ( "hierarchized values = " << values );
//  precicePrint ( "hierarchized oldValues = " << oldValues );
//  precicePrint ( "hierarchized _residual = " << _residual );
//  precicePrint ( "hierarchized residualDelta = " << residualDelta );

  // Compute and perform relaxation with aitken factor
  nominators[0] = _residual(0) * residualDelta(0) +
                  _residual(entries-1) * residualDelta(entries-1);
  denominators[0] = residualDelta(0) * residualDelta(0) +
                    residualDelta(entries-1) * residualDelta(entries-1);
  computeAitkenFactor ( 0, nominators[0], denominators[0] );
  double omega = _aitkenFactors[0];
  double oneMinusOmega = 1.0 - omega;
  for (DataMap::value_type &pair : cplData) {
    auto& values = *pair.second->values;
    const auto& oldValues = pair.second->oldValues.col(0);
    values(0) = values(0) * omega + oldValues(0) * oneMinusOmega;
    values(entries-1) = values(entries-1) * omega + oldValues(entries-1) * oneMinusOmega;
  }
  treatedEntries = 2;
  entriesCurrentLevel = 1;
  for ( size_t level=1; level < _aitkenFactors.size(); level++ ) {
    size_t stepsize = (entries - 1) / std::pow(2.0, (int)(level-1));
    size_t index = stepsize / 2;
    for ( size_t i=0; i < entriesCurrentLevel; i++ ) {
      nominators[level] += _residual(index) * residualDelta(index);
      denominators[level] += residualDelta(index) * residualDelta(index);
      index += stepsize;
    }
    computeAitkenFactor ( level, nominators[level], denominators[level] );
    omega = _aitkenFactors[level];
    oneMinusOmega = 1.0 - omega;
    //for ( DataMap::value_type & pair : cplData ) {
    //  DataValues & values = *pair.second.values;
    //  DataValues & oldValues = pair.second.oldValues.getColumn(0);
      index = stepsize / 2;
      for ( size_t i=0; i < entriesCurrentLevel; i++ ) {
        values(index) = values(index) * omega + oldValues(index) * oneMinusOmega;
        index += stepsize;
      }
    //}
    treatedEntries += entriesCurrentLevel;
    entriesCurrentLevel *= 2;
  }
  assertion ( treatedEntries == entries );

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

//  precicePrint ( "relaxed hierarchized values = " << values );


  // Dehierarchization
  treatedEntries = 2;
  entriesCurrentLevel = 1;
  for ( size_t level=1; level < _aitkenFactors.size(); level++ ) {
    size_t stepsize = (entries - 1) / std::pow(2.0, (int)(level-1));
    size_t index = stepsize / 2;
    for ( size_t i=0; i < entriesCurrentLevel; i++ ) {
      values(index) +=
          (values(index - stepsize/2) + values(index + stepsize/2)) / 2.0;
      oldValues(index) +=
          (oldValues(index - stepsize/2) + oldValues(index + stepsize/2)) / 2.0;
      index += stepsize;
    }
    treatedEntries += entriesCurrentLevel;
    entriesCurrentLevel *= 2;
  }
  assertion ( treatedEntries == entries );
//  precicePrint ( "relaxed values = " << values );
//  precicePrint ( "oldValues = " << oldValues );


  // save back oldValues in cplData. Eigen does not allow call by ref for blocks (cols, rows), thus explicitly write back
  cplData[*_dataIDs.begin()]->oldValues.col(0) = oldValues;

  _iterationCounter ++;
}

//void HierarchicalAitkenPostProcessing:: performPostProcessing
//(
//  DataMap & cplData )
//{
//  preciceTrace ( "performPostProcessing()" );
//  typedef utils::DynVector DataValues;
//
//  // Compute aitken relaxation factor
//  assertion ( utils::contained(_dataID, cplData) );
//  DataValues & values = *cplData[_dataID].values;
//  DataValues & oldValues = cplData[_dataID].oldValues.getColumn(0);
//
//  // Compute current residuals
//  DataValues residual ( values );
//  residual -= oldValues;
//
//  // Compute residual deltas and temporarily store it in _residuals
//  DataValues residualDelta ( _residual );
//  residualDelta *= -1.0;
//  residualDelta += residual;
//
//  std::vector<double> nominators ( _aitkenFactors.size(), 0.0 );
//  std::vector<double> denominators ( _aitkenFactors.size(), 0.0 );
//
//  precicePrint ( "values = " << values );
//  precicePrint ( "oldValues = " << oldValues );
//  precicePrint ( "_residual = " << _residual );
//  precicePrint ( "residualDelta = " << residualDelta );
//
//  // Hierarchize entries and compute aitken factors
//  size_t entries = residual.size ();
//  size_t indexCenter = (entries+1) / 2;
//  nominators[0] = _residual[indexCenter] * residualDelta[indexCenter];
//  denominators[0] = residualDelta[indexCenter] * residualDelta[indexCenter];
//  computeAitkenFactor ( 0, nominators[0], denominators[0] );
//  size_t treatedEntries = 1;
//  size_t entriesCurrentLevel = std::pow(2.0, _aitkenFactors.size() - 1);
//  for ( size_t level=_aitkenFactors.size()-1; level > 0; level-- ) {
//    size_t stepsize = (entries + 1) / std::pow(2.0, level);
//    precicePrint ( "Level = " << level << ", treatedEntries = " << treatedEntries
//                   << ", stepsize = " << stepsize );
//    assertion ( stepsize % 2 == 0 );
//
//    // Leftmost entry
//    size_t indexLeft = stepsize/2 - 1;
//    _residual[indexLeft] -= _residual[indexLeft + stepsize/2] / 2.0;
//    residualDelta[indexLeft] -= residualDelta[indexLeft + stepsize/2] / 2.0;
//    values[indexLeft] -= values[indexLeft + stepsize/2] / 2.0;
//    oldValues[indexLeft] -= oldValues[indexLeft + stepsize/2] / 2.0;
//    nominators[level] += _residual[indexLeft] * residualDelta[indexLeft];
//    denominators[level] += residualDelta[indexLeft] * residualDelta[indexLeft];
//
//    // Middle entries
//    size_t index = indexLeft + stepsize;
//    for ( size_t i=2; i < entriesCurrentLevel; i++ ) {
////      precicePrint ( "2" );
//      _residual[index] -= ( _residual[index - stepsize/2] +
//                            _residual[index + stepsize/2] ) / 2.0;
//      residualDelta[index] -= ( residualDelta[index - stepsize/2] +
//                                residualDelta[index + stepsize/2] ) / 2.0;
//      values[index] -= ( values[index - stepsize/2] +
//                         values[index + stepsize/2] ) / 2.0;
//      oldValues[index] -= ( oldValues[index - stepsize/2] +
//                            oldValues[index + stepsize/2] ) / 2.0;
//      nominators[level] += _residual[index] * residualDelta[index];
//      denominators[level] += residualDelta[index] * residualDelta[index];
//      index += stepsize;
//    }
//
//    // Rightmost entry
//    size_t indexRight = entries - stepsize/2;
//    _residual[indexRight] -= _residual[indexRight - stepsize/2] / 2.0;
//    residualDelta[indexRight] -= residualDelta[indexRight - stepsize/2] / 2.0;
//    values[indexRight] -= values[indexRight - stepsize/2] / 2.0;
//    oldValues[indexRight] -= oldValues[indexRight - stepsize/2] / 2.0;
//    nominators[level] += _residual[indexRight] * residualDelta[indexRight];
//    denominators[level] += residualDelta[indexRight] * residualDelta[indexRight];
//
//    computeAitkenFactor ( level, nominators[level], denominators[level] );
//    treatedEntries += entriesCurrentLevel;
//    entriesCurrentLevel /= 2;
//  }
//  assertion ( treatedEntries == entries );
//  precicePrint ( "hierarchized values = " << values );
//  precicePrint ( "hierarchized oldValues = " << oldValues );
//  precicePrint ( "hierarchized _residual = " << _residual );
//  precicePrint ( "hierarchized residualDelta = " << residualDelta );
//
//  _residual = residual;
//
//  // Perform relaxation with aitken factor
//  treatedEntries = 0;
//  entriesCurrentLevel = 1;
//  for ( size_t level=0; level < _aitkenFactors.size(); level++ ) {
//    double omega = _aitkenFactors[level];
//    double oneMinusOmega = 1.0 - omega;
//    for ( DataMap::value_type & pair : cplData ) {
//      DataValues & values = *pair.second.values;
//      DataValues & oldValues = pair.second.oldValues.getColumn(0);
//      size_t stepsize = (entries + 1) / std::pow(2.0, level);
//      size_t index = stepsize / 2 - 1;
//      for ( size_t i=0; i < entriesCurrentLevel; i++ ) {
//        values[index] = values[index] * omega + oldValues[index] * oneMinusOmega;
//        index += stepsize;
//      }
//    }
//    treatedEntries += entriesCurrentLevel;
//    entriesCurrentLevel *= 2;
//  }
//  assertion ( treatedEntries == entries );
//  precicePrint ( "relaxed hierarchized values = " << values );
//
//  // Dehierarchization
//  treatedEntries = 1;
//  entriesCurrentLevel = 2;
//  for ( size_t level=1; level < _aitkenFactors.size(); level++ ) {
//    size_t stepsize = (entries + 1) / std::pow(2.0, level);
//    size_t index = stepsize / 2 - 1;
//    values[index] += values[index + stepsize/2] / 2.0;
//    oldValues[index] += oldValues[index + stepsize/2] / 2.0;
//    for ( size_t i=1; i < entriesCurrentLevel - 1; i++ ) {
//      index += stepsize;
//      values[index] +=
//          (values[index - stepsize/2] + values[index + stepsize/2]) / 2.0;
//      oldValues[index] +=
//          (oldValues[index - stepsize/2] + oldValues[index + stepsize/2]) / 2.0;
//    }
//    index += stepsize;
//    values[index] += values[index - stepsize/2] / 2.0;
//    oldValues[index] += oldValues[index - stepsize/2] / 2.0;
//    treatedEntries += entriesCurrentLevel;
//    entriesCurrentLevel *= 2;
//  }
//  assertion ( treatedEntries == entries );
//  precicePrint ( "relaxed values = " << values );
//  precicePrint ( "oldValues = " << oldValues );
//
//  _iterationCounter ++;
//}

void HierarchicalAitkenPostProcessing:: iterationsConverged
(
  DataMap & cplData )
{
  _iterationCounter = 0;
  _residual = Eigen::VectorXd::Constant(_residual.size(), std::numeric_limits<double>::max ());
}

void HierarchicalAitkenPostProcessing:: computeAitkenFactor
(
  size_t level,
  double nominator,
  double denominator )
{
  using namespace tarch::la;
  // Select/compute aitken factor depending on current iteration count
  if ( _iterationCounter == 0 ) {
    //precicePrint ( "First iteration (nom = " << nominator << ", den = " << denominator );
    _aitkenFactors[level] = sign(_aitkenFactors[level]) * min(
        Vector<2,double>(_initialRelaxation, std::abs(_aitkenFactors[level])));
  }
  else {
    //precicePrint ( "Aitken factor = - " << _aitkenFactors[level] << " * "
    //               << nominator << " / " << denominator );
    if ( equals(std::sqrt(denominator), 0.0) ) {
      _aitkenFactors[level] = 1.0;
    }
    else {
      _aitkenFactors[level] = -_aitkenFactors[level] * (nominator / denominator);
    }
  }
  //precicePrint ( "Level " << level << " relaxation factor = " << _aitkenFactors[level] );
}

/** ---------------------------------------------------------------------------------------------
 *         getDesignSpecification()
 *
 * @brief: Returns the design specification corresponding to the given coupling data.
 *         This information is needed for convergence measurements in the coupling scheme.
 *  ---------------------------------------------------------------------------------------------
 */
std::map<int, Eigen::VectorXd> HierarchicalAitkenPostProcessing::getDesignSpecification
(
  DataMap& cplData)
{
  preciceError(__func__, "design specification for Aitken relaxation is not supported yet.");

  std::map<int, Eigen::VectorXd> designSpecifications;
  int off = 0;
  for (int id : _dataIDs) {
      int size = cplData[id]->values->size();
      Eigen::VectorXd q = Eigen::VectorXd::Zero(size);
      for (int i = 0; i < size; i++) {
        q(i) = _designSpecification(i+off);
      }
      off += size;
      std::map<int, Eigen::VectorXd>::value_type pair = std::make_pair(id, q);
      designSpecifications.insert(pair);
    }
  return designSpecifications;
}

void HierarchicalAitkenPostProcessing::setDesignSpecification(
     Eigen::VectorXd& q)
 {
   _designSpecification = q;
   preciceError(__func__, "design specification for Aitken relaxation is not supported yet.");
 }

}}} // namespace precice, cplscheme, impl

