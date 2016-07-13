#pragma once

#include "Dimensions.hpp"

#include "com/Communication.hpp"

#include "tarch/logging/Log.h"

namespace precice {
namespace utils {

/**
 * @brief Utility class for managing Master-Slave operations.
 */
class MasterSlave
{
public:

  static int _rank;

  /// Number of ranks. This includes ranks from both participants, e.g. minimal size is 2.
  static int _size;

  /// The rank of the master. At this time it is hardcodes to 0.
  static int _masterRank;

  /**
   * @brief True if this process is running the master.
   */
  static bool _masterMode;
  /**
   * @brief True if this process is running a slave.
   */
  static bool _slaveMode;

  /**
   * @brief Communication between the master and all slaves.
   */
  static com::Communication::SharedPointer _communication;


  /**
   * @brief Configure the master-slave communication.
   */
  static void configure(int rank, int size);

  /**
   * @brief the l2 norm of a vector is calculated on distributed data.
   */
  static double l2norm(const DynVector& vec);
  static double l2norm(const EigenVector& vec);

  /**
   * @brief the dot product of 2 vectors is calculated on distributed data.
   */
  static double dot(const DynVector& vec1, const DynVector& vec2);
  static double dot(const EigenVector& vec1, const EigenVector& vec2);

  static void reset();

  static void reduceSum(double* sendData, double* rcvData, int size);

  static void reduceSum(int& sendData, int& rcvData, int size);

  static void allreduceSum(double* sendData, double* rcvData, int size);

  static void allreduceSum(double& sendData, double& rcvData, int size);

  static void allreduceSum(int& sendData, int& rcvData, int size);

  static void broadcast(bool& value);

  static void broadcast(double& value);
  
  static void broadcast(double* values, int size);

private:

  static tarch::logging::Log _log;

};


}} // namespace precice, utils
