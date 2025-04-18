#ifndef PRECICE_NO_MPI

#pragma once

#include <mpi.h>
#include <string>

#include "com/Communication.hpp"
#include "com/SharedPointer.hpp"
#include "logging/Logger.hpp"
#include "precice/impl/Types.hpp"

namespace precice::com {
/**
 * @brief Provides implementation for basic MPI point-to-point communication.
 *
 * The methods for establishing a connection between two coupling participants
 * are not implemented and left to subclasses.
 */
class MPICommunication : public ::precice::com::Communication {
public:
  MPICommunication();

  /// Destructor, empty.
  ~MPICommunication() override = default;

  /**
   * @brief Sends a std::string to process with given rank.
   *
   * Default MPI point-to-point communication is used.
   */
  void send(std::string const &itemToSend, Rank rankReceiver) override;

  /// Sends an array of integer values.
  void send(precice::span<const int> itemsToSend, Rank rankReceiver) override;

  /// Asynchronously sends an array of integer values.
  PtrRequest aSend(precice::span<const int> itemsToSend, Rank rankReceiver) override;

  /// Sends an array of double values.
  void send(precice::span<const double> itemsToSend, Rank rankReceiver) override;

  /// Asynchronously sends an array of double values.
  PtrRequest aSend(precice::span<const double> itemsToSend, Rank rankReceiver) override;

  /**
   * @brief Sends a double to process with given rank.
   *
   * Default MPI point-to-point communication is used.
   */
  void send(double itemToSend, Rank rankReceiver) override;

  /// Asynchronously sends a double to process with given rank.
  PtrRequest aSend(const double &itemToSend, Rank rankReceiver) override;

  /**
   * @brief Sends an int to process with given rank.
   *
   * Default MPI point-to-point communication is used.
   */
  void send(int itemToSend, Rank rankReceiver) override;

  /// Asynchronously sends an int to process with given rank.
  PtrRequest aSend(const int &itemToSend, Rank rankReceiver) override;

  /**
   * @brief Sends a bool to process with given rank.
   *
   * Default MPI point-to-point communication is used.
   */
  void send(bool itemToSend, Rank rankReceiver) override;

  /// Asynchronously sends a bool to process with given rank.
  PtrRequest aSend(const bool &itemToSend, Rank rankReceiver) override;

  /**
   * @brief Receives a std::string from process with given rank.
   *
   * Default MPI point-to-point communication is used.
   */
  void receive(std::string &itemToReceive, Rank rankSender) override;

  /// Receives an array of integer values.
  void receive(precice::span<int> itemsToReceive, Rank rankSender) override;

  /// Receives an array of double values.
  void receive(precice::span<double> itemsToReceive, Rank rankSender) override;

  /// Asynchronously receives an array of double values.
  PtrRequest aReceive(precice::span<double> itemsToReceive, int rankSender) override;

  /**
   * @brief Receives a double from process with given rank.
   *
   * Default MPI point-to-point communication is used.
   */
  void receive(double &itemToReceive, Rank rankSender) override;

  /// Asynchronously receives a double from process with given rank.
  PtrRequest aReceive(double &itemToReceive, Rank rankSender) override;

  /**
   * @brief Receives an int from process with given rank.
   *
   * Default MPI point-to-point communication is used.
   */
  void receive(int &itemToReceive, Rank rankSender) override;

  /// Asynchronously receives an int from process with given rank.
  PtrRequest aReceive(int &itemToReceive, Rank rankSender) override;

  /**
   * @brief Receives a bool from process with given rank.
   *
   * Default MPI point-to-point communication is used.
   */
  void receive(bool &itemToReceive, Rank rankSender) override;

  /// Asynchronously receives a bool from process with given rank.
  PtrRequest aReceive(bool &itemToReceive, Rank rankSender) override;

protected:
  /// Returns the communicator.
  virtual MPI_Comm &communicator(Rank rank) = 0;

  virtual Rank rank(int rank) = 0;

private:
  logging::Logger _log{"com::MPICommunication"};
};
} // namespace precice::com

#endif // not PRECICE_NO_MPI
