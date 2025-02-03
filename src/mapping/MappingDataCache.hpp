#pragma once

#include <Eigen/Core>
#include "mapping/SharedPointer.hpp"
#include "mesh/Mesh.hpp"
#include "precice/impl/Types.hpp"
#include "time/Sample.hpp"

namespace precice {
namespace mapping {

/**
 * Struct to store temporary data tied to a specific Mapping-Data pair
 * The usual mapping data structures have no notion about the underlying data.
 * The Cache stores in addition a timestep the data corresponds to.
 */
struct MappingDataCache {
public:
  /// Default constructor
  explicit MappingDataCache(int dataDim);

  // For PUM, we need a sequence of polynomial contributions
  // all data here is stored as a Matrix to encode the number of components
  // the layout is
  // Eigen::MatrixXd (vertices, components);
  // Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor, Eigen::Dynamic, 3> instead of MatrixXd?
  // an alternative layout of this data could be to wrap this in another struct ClusterData which is then within a vector
  std::vector<Eigen::MatrixXd> polynomialContributions;
  // ...and a vector of P/lambdas
  // For conservative mappings, this vector refers to the Au, potentially as a vector components
  std::vector<Eigen::MatrixXd> p;

  // preprocessed data in \ref polynomialContributions or \ref p
  // std::unique_ptr<time::Sample> inData;
  // @todo: We could make this a sample rather than a VectorXd, as we might want to operate on gradient data as well
  // however, there doesn't seem to be any interface to sample data and receiving a sample
  Eigen::VectorXd inData;

  int getDataDimensions() const;

  bool hasDataAtTimeStamp(double time) const;

  // Set the timestamp of the MappingDataCache to the sepcifid time
  void setTimeStamp(double time);

  // Reset all data to zero
  void resetData();

private:
  double    _timeStamp = -1;
  const int _dataDim;
};

// ------------------------------------------------------ HEADER IMPLEMENTATION
inline MappingDataCache::MappingDataCache(int dataDim)
    : _dataDim(dataDim)
{
}

inline int MappingDataCache::getDataDimensions() const
{
  return _dataDim;
}

inline bool MappingDataCache::hasDataAtTimeStamp(double time) const
{
  return math::equals(_timeStamp, time);
}

inline void MappingDataCache::setTimeStamp(double time)
{
  _timeStamp = time;
}

inline void MappingDataCache::resetData()
{
  inData.setZero();

  for (auto &mat : polynomialContributions) {
    mat.setZero();
  }
  for (auto &mat : p) {
    mat.setZero();
  }
}
} // namespace mapping
} // namespace precice
