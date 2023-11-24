#pragma once

#include <Eigen/Core>
#include <map>
#include <vector>
#include "acceleration/Acceleration.hpp"
#include "acceleration/BaseQNAcceleration.hpp"
#include "acceleration/impl/SharedPointer.hpp"

namespace precice {
namespace acceleration {

/**
 * @brief Interface quasi-Newton with interface least-squares approximation.
 *
 * Performs a quasi-Newton to accelerate the convergence of implicit coupling
 * iterations. A least-squares approximation is used to find the best linear
 * combination of old data values. After every coupling iteration, the data
 * values used are enhanced by the new coupling iterates.
 *
 * If more coupling data is present than used to compute the IQN acceleration,
 * this data is relaxed using the same linear combination as computed for the
 * IQN-related data. The data is called "secondary" henceforth and additional
 * old value and data matrices are needed for it.
 */
class IQNILSAcceleration : public BaseQNAcceleration {
public:
  IQNILSAcceleration(
      double                  initialRelaxation,
      bool                    forceInitialRelaxation,
      int                     maxIterationsUsed,
      int                     pastTimeWindowsReused,
      int                     filter,
      double                  singularityLimit,
      std::vector<int>        dataIDs,
      impl::PtrPreconditioner preconditioner);

  virtual ~IQNILSAcceleration() {}

  /// Initializes the acceleration.
  virtual void initialize(const DataMap &cplData);

  /**
    * @brief Marks a iteration sequence as converged.
    *
    * called by the iterationsConverged() method in the BaseQNAcceleration class
    * handles the acceleration specific action after the convergence of one iteration
    */
  virtual void specializedIterationsConverged(const DataMap &cplData);

private:
  /// Secondary data solver output from last iteration.
  std::map<int, precice::time::Storage> _secondaryOldXTildesW;

  // @brief Secondary data x-tilde deltas.
  //
  // Stores x-tilde deltas for data not involved in least-squares computation.

  std::map<int, std::vector<precice::time::Storage>> _secondaryWaveformW;
  std::map<int, std::vector<precice::time::Storage>> _secondaryWaveformWBackup;

  /// updates the secondary W waveforms
  void addSecondaryWaveforms(const DataMap &cplData);

  /// updates the V, W matrices (as well as the matrices for the secondary data)
  virtual void updateDifferenceMatrices(const DataMap &cplData);

  /// computes the IQN-ILS update using QR decomposition
  virtual void computeQNUpdate(const DataMap &cplData);

  /// computes underrelaxation for the secondary data
  virtual void computeUnderrelaxationSecondaryData(const DataMap &cplData);

  /// Removes one iteration from V,W matrices and adapts _matrixCols.
  virtual void removeMatrixColumn(int columnIndex);
};
} // namespace acceleration
} // namespace precice
