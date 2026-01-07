#pragma once

#include <Eigen/Core>
#include <map>
#include <vector>

#include "cplscheme/BaseCouplingScheme.hpp"
#include "cplscheme/SharedPointer.hpp"

namespace precice::io {
class TXTWriter;
class TXTReader;
class TXTTableWriter;
} // namespace precice::io

namespace precice::acceleration {

class Acceleration {
public:
  static const int NOFILTER      = 0;
  static const int QR1FILTER     = 1;
  static const int QR1FILTER_ABS = 2;
  static const int QR2FILTER     = 3;
  static const int PODFILTER     = 4;
  static const int QR3FILTER     = 5;

  /// Map from data ID to data values.
  using DataMap = std::map<int, cplscheme::PtrCouplingData>;

  virtual ~Acceleration() = default;

  virtual std::vector<int> getPrimaryDataIDs() const = 0;

  virtual void initialize(const DataMap &cpldata) = 0;

  virtual void performAcceleration(DataMap &cpldata, double windowStart, double windowEnd) = 0;

  virtual void iterationsConverged(const DataMap &cpldata, double windowStart) = 0;

  virtual void exportState(io::TXTWriter &writer) {}

  virtual void importState(io::TXTReader &reader) {}

  /// Allows the acceleration to add custom columns to the iteration log file
  virtual void addLogEntries(io::TXTTableWriter &writer) const
  {
    // Default: no extra iteration columns
  }

  /// Allows the acceleration to write values to custom iteration log columns
  virtual void writeLogEntries(io::TXTTableWriter &writer) const
  {
    // Default: nothing to write
  }

protected:
  /// Checks if all dataIDs are contained in cplData
  void checkDataIDs(const DataMap &cplData) const;

  /// performs a relaxation given a relaxation factor omega
  static void applyRelaxation(double omega, DataMap &cplData, double windowStart);
};
} // namespace precice::acceleration
