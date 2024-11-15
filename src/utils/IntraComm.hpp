#pragma once

#include <Eigen/Core>

#include "boost/range/irange.hpp"
#include "com/SharedPointer.hpp"
#include "logging/Logger.hpp"
#include "precice/impl/Types.hpp"
#include "precice/span.hpp"

namespace precice {
namespace logging {
class Logger;
} // namespace logging

namespace utils {

/// Utility class for managing intra-participant communication operations.
class IntraComm {
public:
  /// Configures the intra-participant communication.
  static void configure(Rank rank, int size);

  /// Current rank
  static Rank getRank();

  /// Number of ranks. This includes ranks from both participants, e.g. minimal size is 2.
  static int getSize();

  /// Intra-participant communication.
  static com::PtrCommunication &getCommunication()
  {
    return _communication;
  }

  /// Returns an iterable range over salve ranks [1, _size)
  static auto allSecondaryRanks()
  {
    return boost::irange(1, _size);
  }

  /// Returns an iterable range over all ranks [0, _size)
  static auto allRanks()
  {
    return boost::irange(0, _size);
  }

  /// True if this process is running the primary rank.
  static bool isPrimary();

  /// True if this process is running a secondary rank.
  static bool isSecondary();

  /// True if this process is running in parallel
  static bool isParallel();

  /// The l2 norm of a vector is calculated on distributed data.
  static double l2norm(const Eigen::VectorXd &vec);

  // The dot product of 2 vectors is calculated on distributed data.
  static double dot(const Eigen::VectorXd &vec1, const Eigen::VectorXd &vec2);

  static void reset();

  static void reduceSum(precice::span<const double> sendData, precice::span<double> rcvData);

  static void reduceSum(const int &sendData, int &rcvData);

  static void reduceSum(const double &sendData, double &rcvData);

  static void allreduceSum(precice::span<const double> sendData, precice::span<double> rcvData);

  static void allreduceSum(double &sendData, double &rcvData);

  static void allreduceSum(int &sendData, int &rcvData);

  static void broadcast(bool &value);

  static void broadcast(double &value);

  static void broadcast(int &value);

  static void broadcast(precice::span<double> values);

  /** Synchronizes all ranks if syncMode is enabled
   * @see precice::syncMode
   * @see barrier()
   */
  static void synchronize();

  /** Does @synchronize have an effect?
   * @see precice::syncMode
   */
  static bool willSynchronize();

  /// Synchronizes all ranks
  static void barrier();

private:
  static logging::Logger _log;

  /// Current rank
  static Rank _rank;

  /// Number of ranks. This includes ranks from both participants, e.g. minimal size is 2.
  static int _size;

  /// True if this process is running the primary rank.
  static bool _isPrimaryRank;

  /// True if this process is running a secondary rank.
  static bool _isSecondaryRank;

  /// Intra-participant communication.
  static com::PtrCommunication _communication;
};

} // namespace utils
} // namespace precice
