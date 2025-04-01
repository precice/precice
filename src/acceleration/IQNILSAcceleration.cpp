#include "acceleration/IQNILSAcceleration.hpp"
#include <Eigen/Core>
#include <algorithm>
#include <deque>
#include <memory>
#include <utility>

#include "acceleration/impl/Preconditioner.hpp"
#include "acceleration/impl/QRFactorization.hpp"
#include "acceleration/impl/SharedPointer.hpp"
#include "com/Communication.hpp"
#include "com/SharedPointer.hpp"
#include "cplscheme/CouplingData.hpp"
#include "cplscheme/SharedPointer.hpp"
#include "logging/LogMacros.hpp"
#include "utils/EigenHelperFunctions.hpp"
#include "utils/Helpers.hpp"
#include "utils/IntraComm.hpp"
#include "utils/assertion.hpp"

// #include "utils/NumericalCompare.hpp"

using precice::cplscheme::PtrCouplingData;

namespace precice::acceleration {

IQNILSAcceleration::IQNILSAcceleration(
    double                  initialRelaxation,
    bool                    forceInitialRelaxation,
    int                     maxIterationsUsed,
    int                     pastTimeWindowsReused,
    int                     filter,
    double                  singularityLimit,
    std::vector<int>        dataIDs,
    impl::PtrPreconditioner preconditioner,
    bool                    reducedTimeGrid)
    : BaseQNAcceleration(initialRelaxation, forceInitialRelaxation, maxIterationsUsed, pastTimeWindowsReused,
                         filter, singularityLimit, std::move(dataIDs), std::move(preconditioner), reducedTimeGrid)
{
}

void IQNILSAcceleration::updateDifferenceMatrices(
    const DataMap &cplData)
{
  // call the base method for common update of V, W matrices
  BaseQNAcceleration::updateDifferenceMatrices(cplData);
}

void IQNILSAcceleration::computeQNUpdate(Eigen::VectorXd &xUpdate)
{
  PRECICE_TRACE();
  PRECICE_DEBUG("   Compute Newton factors");

  // Calculate QR decomposition of matrix V and solve Rc = -Qr
  Eigen::VectorXd c;

  // for procs with no vertices,
  // qrV.cols() = getLSSystemCols() and _qrV.rows() = 0
  auto Q = _qrV.matrixQ();
  auto R = _qrV.matrixR();

  if (!_hasNodesOnInterface) {
    PRECICE_ASSERT(_qrV.cols() == getLSSystemCols(), _qrV.cols(), getLSSystemCols());
    PRECICE_ASSERT(_qrV.rows() == 0, _qrV.rows());
    PRECICE_ASSERT(Q.size() == 0, Q.size());
  }

  Eigen::VectorXd _local_b = Eigen::VectorXd::Zero(_qrV.cols());
  Eigen::VectorXd _global_b;

  // need to scale the residual to compensate for the scaling in c = R^-1 * Q^T * P^-1 * residual'
  // it is also possible to apply the inverse scaling weights from the right to the vector c
  _preconditioner->apply(_primaryResiduals);
  _local_b = Q.transpose() * _primaryResiduals;
  _preconditioner->revert(_primaryResiduals);
  _local_b *= -1.0; // = -Qr

  PRECICE_ASSERT(c.size() == 0, c.size());
  // reserve memory for c
  utils::append(c, Eigen::VectorXd(Eigen::VectorXd::Zero(_local_b.size())));

  // compute rhs Q^T*res in parallel
  if (!utils::IntraComm::isParallel()) {
    PRECICE_ASSERT(Q.cols() == getLSSystemCols(), Q.cols(), getLSSystemCols());
    // back substitution
    c = R.triangularView<Eigen::Upper>().solve<Eigen::OnTheLeft>(_local_b);
  } else {
    PRECICE_ASSERT(utils::IntraComm::getCommunication() != nullptr);
    PRECICE_ASSERT(utils::IntraComm::getCommunication()->isConnected());
    if (_hasNodesOnInterface) {
      PRECICE_ASSERT(Q.cols() == getLSSystemCols(), Q.cols(), getLSSystemCols());
    }
    PRECICE_ASSERT(_local_b.size() == getLSSystemCols(), _local_b.size(), getLSSystemCols());

    if (utils::IntraComm::isPrimary()) {
      PRECICE_ASSERT(_global_b.size() == 0, _global_b.size());
    }
    utils::append(_global_b, Eigen::VectorXd(Eigen::VectorXd::Zero(_local_b.size())));

    // do a reduce operation to sum up all the _local_b vectors
    utils::IntraComm::reduceSum(_local_b, _global_b);

    // back substitution R*c = b only on the primary rank
    if (utils::IntraComm::isPrimary()) {
      c = R.triangularView<Eigen::Upper>().solve<Eigen::OnTheLeft>(_global_b);
    }

    // broadcast coefficients c to all secondary ranks
    utils::IntraComm::broadcast(c);
  }

  PRECICE_DEBUG("   Apply Newton factors");
  // compute x updates from W and coefficients c, i.e, xUpdate = c*W
  xUpdate = _matrixW * c;
}

void IQNILSAcceleration::specializedIterationsConverged(
    const DataMap &cplData)
{
  PRECICE_TRACE();
  if (_matrixCols.empty()) {
    PRECICE_WARN("The IQN matrix has no columns.");
  } else {
    if (_matrixCols.front() == 0) { // Did only one iteration
      _matrixCols.pop_front();
    }
  }
}

void IQNILSAcceleration::removeMatrixColumn(
    int columnIndex)
{
  PRECICE_ASSERT(_matrixV.cols() > 1);
  BaseQNAcceleration::removeMatrixColumn(columnIndex);
}
} // namespace precice::acceleration
