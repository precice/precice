#pragma once

#include <Eigen/Core>
#include "mapping/SharedPointer.hpp"
#include "mesh/Mesh.hpp"
#include "precice/impl/Types.hpp"
#include "time/Sample.hpp"

namespace precice {
namespace mapping {
namespace impl {

/**
 * Struct to store temporary data structures tied to a specific Mapping-Data pair.
 * The usual implementations from 'Mapping' have no notion about the data they are
 * operating on. One mapping might even be responsible for multiple read or write
 * data, which are sequentially processed one after another in the ParticipantImpl.
 *
 * The Cache here was added in the context of just-in-time data mappings. Fundamentally,
 * a just-in-time mapping can also be responsible for multiple read and write data.
 * However, the just-in-time mapping was specifically designed for particle-mesh coupling,
 * where a particle solver provides coordinates just-in-time for the mapping. If now
 * each particle solver calls a just-in-time read or write data call, then we would need
 * to perform the full time interpolation and mapping evaluation for each API function call.
 * The cache here allows to store temporary data we can reuse across multiple just-in-time
 * API calls to preserve computational efficiency.
 * In the most basic configuration scenarios (NN mapping), the cache simply holds the latest
 * time interpolant to prevent repeated time interpolation for each function call. For
 * more advanced mappings (PUM), the mapping stores intermediate basis-function coefficients
 * of the resulting from the last time interpolant. The cache data always holds time-specific
 * data for all data sites.
 *
 * Note that the MappingDataCache mostly provides the memory and data structures. All data structures
 * are filled within the configured mapping class.
 * Particularly relevant are
 * - Mapping::initializeMappingDataCache (allocates memory and resizes data)
 * - Mapping::updateMappingDataCache (updates the cache data to a new time sample)
 * - Mapping::completeJustInTimeMapping (for conservative mappings: processes the buffered data)
 * - Mapping::mapConsistentAt (uses the cache to map consistent)
 * - Mapping::mapConservativeAt (uses the cache to buffer conservative data)
 * To associate the data the cache holds to a specific timestamp, the cache also stores in
 * such a timestamp.
 */
struct MappingDataCache {
public:
  /// Constructor
  explicit MappingDataCache(int dataDim);

  // The std::vector<> member variables here are specific to PUM:
  // we store cache data for each cluster:
  // We need a sequence of polynomial contributions
  // all data here is stored as a Matrix to encode the number of components
  // the layout is
  // Eigen::MatrixXd (vertices, components);
  // Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor, Eigen::Dynamic, 3> instead of MatrixXd?
  // an alternative layout of this data could be to wrap this in another struct ClusterData which is then within a vector
  std::vector<Eigen::MatrixXd> polynomialContributions;
  // ...and a vector of P/lambdas (the RBF coefficients)
  // For conservative mappings, this vector refers to the Au, potentially as vector components
  std::vector<Eigen::MatrixXd> p;

  // The inData member is for basic mappings (such as NN) and simply stores
  // the latest evaluation of the waveform relaxation to prevent repeated time
  // interpolation
  // @note: We could make this a sample rather than a VectorXd, as we might want to operate on gradient data as well
  // std::unique_ptr<time::Sample> inData;
  // however, there doesn't seem to be any compatible interface in time::Sample, so we simply use a vector here
  Eigen::VectorXd inData;

  /// Returns the number of data components
  int getDataDimensions() const;

  /// Check, if the current data vector (in \p inData or the std::vectors) hold data of time \p time
  bool hasDataAtTimeStamp(double time) const;

  /// Set the timestamp of the MappingDataCache to the specified time
  void setTimeStamp(double time);

  /// Reset all data containers to zero
  void resetData();

private:
  double    _timeStamp = -1;
  const int _dataDim{};
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
} // namespace impl
} // namespace mapping
} // namespace precice
