#include "BroydenPostProcessing.hpp"
#include <Eigen/Core>
#include "cplscheme/CouplingData.hpp"

namespace precice
{
namespace cplscheme
{
namespace impl
{

// logging::Logger BroydenPostProcessing::
//      _log("cplscheme::impl::BroydenPostProcessing");

BroydenPostProcessing::BroydenPostProcessing(
    double            initialRelaxation,
    bool              forceInitialRelaxation,
    int               maxIterationsUsed,
    int               timestepsReused,
    int               filter,
    double            singularityLimit,
    std::vector<int>  dataIDs,
    PtrPreconditioner preconditioner)
    : BaseQNPostProcessing(initialRelaxation, forceInitialRelaxation, maxIterationsUsed, timestepsReused,
                           filter, singularityLimit, dataIDs, preconditioner),
      _maxColumns(maxIterationsUsed)
{}

void BroydenPostProcessing::initialize(
    DataMap &cplData)
{
  // do common QN post processing initialization
  BaseQNPostProcessing::initialize(cplData);

  size_t entries = _residuals.size();

  _invJacobian    = Eigen::MatrixXd::Zero(entries, entries);
  _oldInvJacobian = Eigen::MatrixXd::Zero(entries, entries);
}

void BroydenPostProcessing::computeUnderrelaxationSecondaryData(
    DataMap &cplData)
{
  // Perform underrelaxation with initial relaxation factor for secondary data
  for (int id : _secondaryDataIDs) {
    PtrCouplingData  data   = cplData[id];
    Eigen::VectorXd &values = *(data->values);
    values *= _initialRelaxation; // new * omg
    Eigen::VectorXd &secResiduals = _secondaryResiduals[id];
    secResiduals                  = data->oldValues.col(0); // old
    secResiduals *= 1.0 - _initialRelaxation;               // (1-omg) * old
    values += secResiduals;                                 // (1-omg) * old + new * omg
  }
}

void BroydenPostProcessing::updateDifferenceMatrices(
    DataMap &cplData)
{
  if (not _firstIteration) {
      _currentColumns++;
  }

  // call the base method for common update of V, W matrices
  BaseQNPostProcessing::updateDifferenceMatrices(cplData);
}

void BroydenPostProcessing::computeQNUpdate(PostProcessing::DataMap &cplData, Eigen::VectorXd &xUpdate)
{
  TRACE();

  DEBUG("currentColumns=" << _currentColumns);
  if (_currentColumns > 1) {
    ERROR("truncated IMVJ no longer supported, needs to be parallelized and datastructures need to be changed to Eigen datastructures.");
    DEBUG("compute update with QR-dec");
    //computeNewtonFactorsQRDecomposition(cplData, xUpdate);
  } else {
    DEBUG("compute update with Broyden");
    // ------------- update inverse Jacobian -----------
    // ------------- Broyden Update
    //
    // J_inv = J_inv_n + (w- J_inv_n*v)*v^T/|v|_l2
    // ----------------------------------------- -------

    Eigen::VectorXd v       = _matrixV.col(0);
    Eigen::VectorXd w       = _matrixW.col(0);
    Eigen::MatrixXd JUpdate = Eigen::MatrixXd::Zero(_invJacobian.rows(), _invJacobian.cols());

    DEBUG("took latest column of V,W");

    double          dotproductV = v.dot(v);
    Eigen::VectorXd tmp         = _oldInvJacobian * v; // J_inv*v
    tmp                         = w - tmp;             // (w-J_inv*v)
    tmp                         = tmp / dotproductV;   // (w-J_inv*v)/|v|_l2
    DEBUG("did step (W-J_inv*v)/|v|");

    assertion(tmp.size() == v.size(), tmp.size(), v.size());
    JUpdate = tmp * v.transpose();
    DEBUG("multiplied (w-J_inv*v)/|v| * v^T");

    assertion(_invJacobian.rows() == JUpdate.rows(), _invJacobian.rows(), JUpdate.rows());
    _invJacobian = _oldInvJacobian + JUpdate;

    // solve delta_x = - J_inv*residuals
    xUpdate = _invJacobian * (-_residuals);
  }

  if (_currentColumns >= _maxColumns) {
    _currentColumns = 0;
    _oldInvJacobian = _invJacobian;
  }
}

void BroydenPostProcessing::specializedIterationsConverged(
    DataMap &cplData)
{
  _currentColumns = 0;
  // store old Jacobian
  _oldInvJacobian = _invJacobian;
}
}
}
} // namespace precice, cplscheme, impl
