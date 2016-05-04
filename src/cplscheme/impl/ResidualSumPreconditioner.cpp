// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "ResidualSumPreconditioner.hpp"
#include "utils/MasterSlave.hpp"

namespace precice {
namespace cplscheme {
namespace impl {

tarch::logging::Log ResidualSumPreconditioner::
   _log ( "precice::cplscheme::ResidualSumPreconditioner" );

ResidualSumPreconditioner:: ResidualSumPreconditioner
(
    std::vector<int> dimensions,
    int maxNonConstTimesteps)
 :
    Preconditioner (dimensions,
         maxNonConstTimesteps)
{
  _residualSum.resize(dimensions.size(),0.0);
}


void ResidualSumPreconditioner::_update_(bool timestepComplete, const Eigen::VectorXd& oldValues, const Eigen::VectorXd& res)
{
  if(not timestepComplete){
    std::vector<double> norms(_dimensions.size(),0.0);

    double sum = 0.0;

    int offset = 0;
    for(size_t k=0; k<_dimensions.size(); k++){
      Eigen::VectorXd part = Eigen::VectorXd::Zero(_dimensions[k]*_sizeOfSubVector);
      for(int i=0; i<_dimensions[k]*_sizeOfSubVector; i++){
        part(i) = res(i+offset);
      }
      norms[k] = utils::MasterSlave::dot(part,part);
      sum += norms[k];
      offset += _dimensions[k]*_sizeOfSubVector;
      norms[k] = std::sqrt(norms[k]);
    }
    sum = std::sqrt(sum);
    assertion(sum>0);

    for(size_t k=0; k<_dimensions.size(); k++){
      _residualSum[k] += norms[k] / sum;
      assertion(_residualSum[k]>0);
    }

    offset = 0;
    for(size_t k=0; k<_dimensions.size(); k++){
      for(int i=0; i<_dimensions[k]*_sizeOfSubVector; i++){
        _weights[i+offset] = 1 / _residualSum[k];
        _invWeights[i+offset] = _residualSum[k];
      }
      offset += _dimensions[k]*_sizeOfSubVector;
    }

    _requireNewQR = true;
  }
  else{
    for(size_t k=0; k<_dimensions.size(); k++){
      _residualSum[k] = 0.0;
    }
  }
}

}}} // namespace precice, cplscheme
