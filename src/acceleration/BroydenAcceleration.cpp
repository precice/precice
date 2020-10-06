#include "acceleration/BroydenAcceleration.hpp"
#include <Eigen/Core>
#include <algorithm>
#include <map>
#include <memory>
#include <ostream>
#include <stddef.h>
#include "acceleration/impl/QRFactorization.hpp"
#include "cplscheme/CouplingData.hpp"
#include "cplscheme/SharedPointer.hpp"
#include "logging/LogMacros.hpp"
#include "utils/assertion.hpp"

namespace precice {
namespace acceleration {

using namespace precice::acceleration::impl;

BroydenAcceleration::BroydenAcceleration(
    double                  initialRelaxation,
    bool                    forceInitialRelaxation,
    int                     maxIterationsUsed,
    int                     timestepsReused,
    int                     filter,
    double                  singularityLimit,
    std::vector<int>        dataIDs,
    impl::PtrPreconditioner preconditioner)
    : BaseQNAcceleration(initialRelaxation, forceInitialRelaxation, maxIterationsUsed, timestepsReused,
                         filter, singularityLimit, dataIDs, preconditioner),
      _maxColumns(maxIterationsUsed)
{
}

void BroydenAcceleration::initialize(
    DataMap &cplData)
{
  // do common QN acceleration initialization
  BaseQNAcceleration::initialize(cplData);

  size_t entries = _residuals.size();

  _invJacobian    = Eigen::MatrixXd::Zero(entries, entries);
  _oldInvJacobian = Eigen::MatrixXd::Zero(entries, entries);
}

void BroydenAcceleration::computeUnderrelaxationSecondaryData(
    DataMap &cplData)
{
  // Perform underrelaxation with initial relaxation factor for secondary data
  for (int id : _secondaryDataIDs) {
    cplscheme::PtrCouplingData data   = cplData[id];
    Eigen::VectorXd &          values = data->values();
    values *= _initialRelaxation; // new * omg
    Eigen::VectorXd &secResiduals = _secondaryResiduals[id];
    secResiduals                  = data->oldValues.col(0); // old
    secResiduals *= 1.0 - _initialRelaxation;               // (1-omg) * old
    values += secResiduals;                                 // (1-omg) * old + new * omg
  }
}

void BroydenAcceleration::updateDifferenceMatrices(
    DataMap &cplData)
{
  if (not _firstIteration) {
    _currentColumns++;
  }

  // call the base method for common update of V, W matrices
  BaseQNAcceleration::updateDifferenceMatrices(cplData);
}

void BroydenAcceleration::computeQNUpdate(Acceleration::DataMap &cplData, Eigen::VectorXd &xUpdate)
{
  PRECICE_TRACE();

  PRECICE_DEBUG("currentColumns=" << _currentColumns);
  if (_currentColumns > 1) {
    PRECICE_ERROR("Truncated IMVJ is no longer supported. Please use IMVJ with restart mode instead.");
    PRECICE_DEBUG("compute update with QR-dec");
    //computeNewtonFactorsQRDecomposition(cplData, xUpdate);
  } else {
    PRECICE_DEBUG("compute update with Broyden");
    // ------------- update inverse Jacobian -----------
    // ------------- Broyden Update
    //
    // J_inv = J_inv_n + (w- J_inv_n*v)*v^T/|v|_l2
    // ----------------------------------------- -------

    Eigen::VectorXd v       = _matrixV.col(0);
    Eigen::VectorXd w       = _matrixW.col(0);
    Eigen::MatrixXd JUpdate = Eigen::MatrixXd::Zero(_invJacobian.rows(), _invJacobian.cols());

    PRECICE_DEBUG("took latest column of V,W");

    double          dotproductV = v.dot(v);
    Eigen::VectorXd tmp         = _oldInvJacobian * v; // J_inv*v
    tmp                         = w - tmp;             // (w-J_inv*v)
    tmp                         = tmp / dotproductV;   // (w-J_inv*v)/|v|_l2
    PRECICE_DEBUG("did step (W-J_inv*v)/|v|");

    PRECICE_ASSERT(tmp.size() == v.size(), tmp.size(), v.size());
    JUpdate = tmp * v.transpose();
    PRECICE_DEBUG("multiplied (w-J_inv*v)/|v| * v^T");

    PRECICE_ASSERT(_invJacobian.rows() == JUpdate.rows(), _invJacobian.rows(), JUpdate.rows());
    _invJacobian = _oldInvJacobian + JUpdate;

    // solve delta_x = - J_inv*residuals
    xUpdate = _invJacobian * (-_residuals);
  }

  if (_currentColumns >= _maxColumns) {
    _currentColumns = 0;
    _oldInvJacobian = _invJacobian;
  }
}

void BroydenAcceleration::specializedIterationsConverged(
    DataMap &cplData)
{
  _currentColumns = 0;
  // store old Jacobian
  _oldInvJacobian = _invJacobian;
}
} // namespace acceleration
} // namespace precice
