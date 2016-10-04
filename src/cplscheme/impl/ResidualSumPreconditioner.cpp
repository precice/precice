#include "ResidualSumPreconditioner.hpp"
#include "utils/MasterSlave.hpp"

namespace precice {
namespace cplscheme {
namespace impl {

logging::Logger ResidualSumPreconditioner::
   _log ( "precice::cplscheme::ResidualSumPreconditioner" );

ResidualSumPreconditioner:: ResidualSumPreconditioner
(
    int maxNonConstTimesteps)
 :
    Preconditioner (
         maxNonConstTimesteps)
{}

void ResidualSumPreconditioner::initialize(std::vector<size_t>& svs){
  preciceTrace("initialize()");
  Preconditioner::initialize(svs);

  _residualSum.resize(_subVectorSizes.size(),0.0);
}


void ResidualSumPreconditioner::_update_(bool timestepComplete, const Eigen::VectorXd& oldValues, const Eigen::VectorXd& res)
{
  if(not timestepComplete){
    std::vector<double> norms(_subVectorSizes.size(),0.0);

    double sum = 0.0;

    int offset = 0;
    for(size_t k=0; k<_subVectorSizes.size(); k++){
      Eigen::VectorXd part = Eigen::VectorXd::Zero(_subVectorSizes[k]);
      for(size_t i=0; i<_subVectorSizes[k]; i++){
        part(i) = res(i+offset);
      }
      norms[k] = utils::MasterSlave::dot(part,part);
      sum += norms[k];
      offset += _subVectorSizes[k];
      norms[k] = std::sqrt(norms[k]);
    }
    sum = std::sqrt(sum);
    assertion(sum>0);

    for(size_t k=0; k<_subVectorSizes.size(); k++){
      _residualSum[k] += norms[k] / sum;
      assertion(_residualSum[k]>0);
    }

    offset = 0;
    for(size_t k=0; k<_subVectorSizes.size(); k++){
      for(size_t i=0; i<_subVectorSizes[k]; i++){
        _weights[i+offset] = 1 / _residualSum[k];
        _invWeights[i+offset] = _residualSum[k];
      }
      offset += _subVectorSizes[k];
    }

    _requireNewQR = true;
  }
  else{
    for(size_t k=0; k<_subVectorSizes.size(); k++){
      _residualSum[k] = 0.0;
    }
  }
}

}}} // namespace precice, cplscheme
