// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "AitkenPostProcessing.hpp"
#include "../CouplingData.hpp"
#include "mesh/Data.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/Mesh.hpp"
#include "tarch/la/DynamicVector.h"
#include "tarch/la/MatrixVectorOperations.h"
#include "utils/Globals.hpp"
#include "utils/Dimensions.hpp"
#include <limits>

namespace precice {
namespace cplscheme {
namespace impl {

tarch::logging::Log AitkenPostProcessing::
  _log ( "precice::cplscheme::AitkenPostProcessing" );

AitkenPostProcessing:: AitkenPostProcessing
(
  double initialRelaxation,
  int    dataID )
:
  PostProcessing (),
  _initialRelaxation ( initialRelaxation ),
  _dataID ( dataID ),
  _aitkenFactor ( initialRelaxation ),
  _iterationCounter ( 0 ),
  _residuals ()
{
  preciceCheck ( (_initialRelaxation > 0.0) && (_initialRelaxation <= 1.0),
                 "AitkenPostProcessing()",
                 "Initial relaxation factor for aitken post processing has to "
                 << "be larger than zero and smaller or equal than one!" );
}

void AitkenPostProcessing:: initialize
(
  DataMap& cpldata )
{
  preciceCheck(utils::contained(_dataID, cpldata), "initialize()",
               "Data with ID " << _dataID
               << " is not contained in data given at initialization!" );
  size_t entries = cpldata[_dataID].values->size();
  assertion(entries > 0);
  double initializer = std::numeric_limits<double>::max();
  utils::DynVector toAppend(entries, initializer);
  _residuals.append(toAppend);

  // Append column for old values if not done by coupling scheme yet
  foreach (DataMap::value_type& pair, cpldata){
    int cols = pair.second.oldValues.cols();
    if (cols < 1){
      assertion1(pair.second.values->size() > 0, pair.first);
      pair.second.oldValues.append(CouplingData::DataMatrix(
        pair.second.values->size(), 1, 0.0));
    }
  }
}

void AitkenPostProcessing:: performPostProcessing
(
  DataMap& cplData )
{
  preciceTrace("performPostProcessing()");
  typedef utils::DynVector DataValues;
  using namespace tarch::la;

  // Compute aitken relaxation factor
  assertion(utils::contained(_dataID, cplData));
  DataValues& values = *cplData[_dataID].values;
  DataValues& oldValues = cplData[_dataID].oldValues.column(0);

  // Compute current residuals
  DataValues residuals(values);
  residuals -= oldValues;

  // Compute residual deltas and temporarily store it in _residuals
  DataValues residualDeltas(_residuals);
  residualDeltas *= -1.0;
  residualDeltas += residuals;

  // Select/compute aitken factor depending on current iteration count
  if (_iterationCounter == 0){
    _aitkenFactor = sign(_aitkenFactor) * min(
                    utils::Vector2D(_initialRelaxation, std::abs(_aitkenFactor)));
  }
  else {
    // compute fraction of aitken factor with residuals and residual deltas
    double nominator = dot(_residuals, residualDeltas);
    double denominator = dot(residualDeltas, residualDeltas);
    _aitkenFactor = -_aitkenFactor * (nominator / denominator);
  }

  // Perform relaxation with aitken factor
  double omega = _aitkenFactor;
  double oneMinusOmega = 1.0 - omega;
  foreach ( DataMap::value_type& pair, cplData ) {
    DataValues& values = *pair.second.values;
    DataValues& oldValues = pair.second.oldValues.column(0);
    values *= omega;
    for ( int i=0; i < values.size(); i++ ) {
      values[i] += oldValues[i] * oneMinusOmega;
    }
  }

  // Store residuals for next iteration
  _residuals = residuals;

  _iterationCounter++;
}

void AitkenPostProcessing:: iterationsConverged
(
  DataMap& cplData )
{
  _iterationCounter = 0;
  assign(_residuals) = std::numeric_limits<double>::max();
}

}}} // namespace precice, cplscheme, impl
