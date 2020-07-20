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

  /// Map from data ID to data values.
  using DataMap   = std::map<int, cplscheme::PtrCouplingData>;
  using ValuesMap = std::map<int, Eigen::VectorXd>;

  virtual ~Acceleration() {}

  virtual std::vector<int> getDataIDs() const = 0;

  virtual void initialize(DataMap &cpldata) = 0;

  virtual void performAcceleration(DataMap &cpldata) = 0;

  virtual void iterationsConverged(DataMap &cpldata) = 0;

  virtual void exportState(io::TXTWriter &writer) {}

  virtual void importState(io::TXTReader &reader) {}

  /// Gives the number of QN columns that where filtered out (i.e. deleted) in this time window
  virtual int getDeletedColumns() const
  {
    return 0;
  }

  /// Gives the number of QN columns that went out of scope in this time window
  virtual int getDroppedColumns() const
  {
    return 0;
  }

  /// Gives the number of current QN columns (LS = least squares)
  virtual int getLSSystemCols() const
  {
    return 0;
  }

protected:
  /// Checks if all dataIDs are contained in cplData
  void checkDataIDs(DataMap const &cplData) const;
};
} // namespace acceleration
} // namespace precice
