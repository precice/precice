#ifndef PRECICE_NO_MPI

#pragma once

#include <mpi.h>
#include <string>
#include <vector>

#include "com/Communication.hpp"
#include "com/SharedPointer.hpp"
#include "logging/Logger.hpp"
#include "precice/types.hpp"

namespace precice {
namespace com {
/**
 * @brief Provides implementation for basic MPI point-to-point communication.
 *
 * The methods for establishing a connection between two coupling participants
 * are not implemented and left to subclasses.
 */
class MPICommunication : public Communication {
public:
  MPICommunication();

  /// Destructor, empty.
  virtual ~MPICommunication()
  {
  }

  /**
   * @brief Sends a std::string to process with given rank.
   *
   * Default MPI point-to-point communication is used.
   */
  virtual void send(std::string const &itemToSend, Rank rankReceiver) override;

  /// Sends an array of integer values.
  virtual void send(const int *itemsToSend, int size, Rank rankReceiver) override;

  /// Asynchronously sends an array of integer values.
  virtual PtrRequest aSend(const int *itemsToSend, int size, Rank rankReceiver) override;

  /// Sends an array of double values.
  virtual void send(const double *itemsToSend, int size, Rank rankReceiver) override;

  /// Asynchronously sends an array of double values.
  virtual PtrRequest aSend(const double *itemsToSend, int size, Rank rankReceiver) override;

  virtual PtrRequest aSend(std::vector<double> const &itemsToSend, Rank rankReceiver) override;

  /**
   * @brief Sends a double to process with given rank.
   *
   * Default MPI point-to-point communication is used.
   */
  virtual void send(double itemToSend, Rank rankReceiver) override;

  /// Asynchronously sends a double to process with given rank.
  virtual PtrRequest aSend(const double &itemToSend, Rank rankReceiver) override;

  /**
   * @brief Sends an int to process with given rank.
   *
   * Default MPI point-to-point communication is used.
   */
  virtual void send(int itemToSend, Rank rankReceiver) override;

  /// Asynchronously sends an int to process with given rank.
  virtual PtrRequest aSend(const int &itemToSend, Rank rankReceiver) override;

  /**
   * @brief Sends a bool to process with given rank.
   *
   * Default MPI point-to-point communication is used.
   */
  virtual void send(bool itemToSend, Rank rankReceiver) override;

  /// Asynchronously sends a bool to process with given rank.
  virtual PtrRequest aSend(const bool &itemToSend, Rank rankReceiver) override;

  /**
   * @brief Receives a std::string from process with given rank.
   *
   * Default MPI point-to-point communication is used.
   */
  virtual void receive(std::string &itemToReceive, Rank rankSender) override;

  /// Receives an array of integer values.
  virtual void receive(int *itemsToReceive, int size, Rank rankSender) override;

  /// Receives an array of double values.
  virtual void receive(double *itemsToReceive, int size, Rank rankSender) override;

  /// Asynchronously receives an array of double values.
  virtual PtrRequest aReceive(double *itemsToReceive,
                              int     size,
                              int     rankSender) override;

  virtual PtrRequest aReceive(std::vector<double> &itemsToReceive, Rank rankSender) override;

  /**
   * @brief Receives a double from process with given rank.
   *
   * Default MPI point-to-point communication is used.
   */
  virtual void receive(double &itemToReceive, Rank rankSender) override;

  /// Asynchronously receives a double from process with given rank.
  virtual PtrRequest aReceive(double &itemToReceive, Rank rankSender) override;

  /**
   * @brief Receives an int from process with given rank.
   *
   * Default MPI point-to-point communication is used.
   */
  virtual void receive(int &itemToReceive, Rank rankSender) override;

  /// Asynchronously receives an int from process with given rank.
  virtual PtrRequest aReceive(int &itemToReceive, Rank rankSender) override;

  /**
   * @brief Receives a bool from process with given rank.
   *
   * Default MPI point-to-point communication is used.
   */
  virtual void receive(bool &itemToReceive, Rank rankSender) override;

  /// Asynchronously receives a bool from process with given rank.
  virtual PtrRequest aReceive(bool &itemToReceive, Rank rankSender) override;

  void send(std::vector<int> const &v, Rank rankReceiver) override;
  void receive(std::vector<int> &v, Rank rankSender) override;

  void send(std::vector<double> const &v, Rank rankReceiver) override;
  void receive(std::vector<double> &v, Rank rankSender) override;

protected:
  /// Returns the communicator.
  virtual MPI_Comm &communicator(Rank rank) = 0;

  virtual Rank rank(int rank) = 0;

private:
  logging::Logger _log{"com::MPICommunication"};
};
} // namespace com
} // namespace precice

#endif // not PRECICE_NO_MPI
