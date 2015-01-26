// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "ConstantRelaxationPostProcessing.hpp"
#include "../CouplingData.hpp"
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"
#include "tarch/la/DynamicVector.h"
#include "utils/Globals.hpp"
#include "utils/Dimensions.hpp"
//#include "utils/NumericalCompare.hpp"

namespace precice {
namespace cplscheme {
namespace impl {

tarch::logging::Log ConstantRelaxationPostProcessing::
  _log("precice::cplscheme::ConstantRelaxationPostProcessing");

ConstantRelaxationPostProcessing:: ConstantRelaxationPostProcessing
(
  double relaxation,
  std::vector<int>    dataIDs )
:
  PostProcessing(),
  _relaxation(relaxation),
  _dataIDs(dataIDs)
{
  preciceCheck((relaxation > 0.0) && (relaxation <= 1.0),
               "ConstantRelaxationPostProcessing()",
               "Relaxation factor for constant relaxation post processing "
               << "has to be larger than zero and smaller or equal than one!");
}

void ConstantRelaxationPostProcessing:: initialize
(
  DataMap& cplData )
{
  preciceCheck(utils::contained(*_dataIDs.begin(), cplData), "initialize()",
               "Data with ID " << *_dataIDs.begin()
               << " is not contained in data given at initialization!");

  // Append column for old values if not done by coupling scheme yet
  foreach (DataMap::value_type& pair, cplData){
    int cols = pair.second->oldValues.cols();
    if (cols < 1){
      assertion1(pair.second->values->size() > 0, pair.first);
      pair.second->oldValues.append(CouplingData::DataMatrix(
        pair.second->values->size(), 1, 0.0));
    }
  }
}

void ConstantRelaxationPostProcessing:: performPostProcessing
(
  DataMap& cplData )
{
  preciceTrace("performPostProcessing()");
  double omega = _relaxation;
  double oneMinusOmega = 1.0 - omega;
  foreach (DataMap::value_type & pair, cplData){
    utils::DynVector& values = * pair.second->values;
    utils::DynVector& oldValues = pair.second->oldValues.column(0);
    values *= omega;
    values += oldValues * oneMinusOmega;
    preciceDebug("pp values" << values);
  }
}

}}} // namespace precice, cplscheme, impl
