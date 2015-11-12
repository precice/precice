// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "ConstantPreconditioner.hpp"

namespace precice {
namespace cplscheme {
namespace impl {

tarch::logging::Log ConstantPreconditioner::
   _log ( "precice::cplscheme::ConstantPreconditioner" );

ConstantPreconditioner:: ConstantPreconditioner
(
   std::vector<int> dimensions,
   std::vector<double> factors)
:
   Preconditioner(dimensions),
   _factors(factors)
{}

void ConstantPreconditioner::initialize(int N){
  preciceTrace("initialize()");
  Preconditioner::initialize(N);

  assertion(_factors.size()==_dimensions.size());

  int offset = 0;
  for(size_t k=0; k<_dimensions.size(); k++){
    for(int i=0; i<_dimensions[k]*_sizeOfSubVector; i++){
      _weights[i+offset] = 1.0 / _factors[k];
      _invWeights[i+offset] = _factors[k];
    }
    offset += _dimensions[k]*_sizeOfSubVector;
  }
}

void ConstantPreconditioner::update(bool timestepComplete, const DataValues& oldValues, const DataValues& res)
{
  preciceTrace("update()");
  //nothing to do here
}

void ConstantPreconditioner::update(bool timestepComplete, const Eigen::VectorXd& oldValues, const Eigen::VectorXd& res)
{
  preciceTrace("update()");
  //nothing to do here
}

}}} // namespace precice, cplscheme
