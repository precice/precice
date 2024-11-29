#pragma once

#include <Eigen/Core>
#include <map>
#include <vector>

#include "cplscheme/BaseCouplingScheme.hpp"
#include "cplscheme/SharedPointer.hpp"

namespace precice {
namespace io {
class TXTWriter;
class TXTReader;
} // namespace io
} // namespace precice

namespace precice {
namespace acceleration {

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

  virtual void performAcceleration(DataMap &cpldata) = 0;

  virtual void iterationsConverged(const DataMap &cpldata) = 0;

  virtual void exportState(io::TXTWriter &writer) {}

  virtual void importState(io::TXTReader &reader) {}

protected:
  /// Checks if all dataIDs are contained in cplData
  void checkDataIDs(const DataMap &cplData) const;

  /// performs a relaxation given a relaxation factor omega
  static void applyRelaxation(double omega, DataMap &cplData);
};
} // namespace acceleration
} // namespace precice
