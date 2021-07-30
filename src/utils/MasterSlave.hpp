#pragma once

#include <Eigen/Core>

#include "boost/range/irange.hpp"
#include "com/SharedPointer.hpp"
#include "logging/Logger.hpp"
#include "precice/types.hpp"

namespace precice {
namespace logging {
class Logger;
} // namespace logging

namespace utils {

/// Utility class for managing Master-Slave operations.
class MasterSlave {
public:
  /// Communication between the master and all slaves.
  static com::PtrCommunication _communication;

  /// Configures the master-slave communication.
  static void configure(Rank rank, int size);

  /// Current rank
  static Rank getRank();

  /// Number of ranks. This includes ranks from both participants, e.g. minimal size is 2.
  static int getSize();

  /// Returns an iterable range over salve ranks [1, _size)
  static auto allSlaves()
  {
    return boost::irange(1, _size);
  }

  /// Returns an iterable range over all ranks [0, _size)
  static auto allRanks()
  {
    return boost::irange(0, _size);
  }

  /// True if this process is running the master.
  static bool isMaster();

  /// True if this process is running a slave.
  static bool isSlave();

  /// True if this process is running in parallel
  static bool isParallel();

  /// The l2 norm of a vector is calculated on distributed data.
  static double l2norm(const Eigen::VectorXd &vec);

  // The dot product of 2 vectors is calculated on distributed data.
  static double dot(const Eigen::VectorXd &vec1, const Eigen::VectorXd &vec2);

  static void reset();

  static void reduceSum(double *sendData, double *rcvData, int size);

  static void reduceSum(int &sendData, int &rcvData, int size);

  static void allreduceSum(double *sendData, double *rcvData, int size);

  static void allreduceSum(double &sendData, double &rcvData, int size);

  static void allreduceSum(int &sendData, int &rcvData, int size);

  static void broadcast(bool &value);

  static void broadcast(double &value);

  static void broadcast(double *values, int size);

private:
  static logging::Logger _log;

  /// Current rank
  static Rank _rank;

  /// Number of ranks. This includes ranks from both participants, e.g. minimal size is 2.
  static int _size;

  /// True if this process is running the master.
  static bool _isMaster;

  /// True if this process is running a slave.
  static bool _isSlave;
};

} // namespace utils
} // namespace precice
