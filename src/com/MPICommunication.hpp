#ifndef PRECICE_NO_MPI

#pragma once

#include <mpi.h>
#include <string>
#include <vector>
#include "com/Communication.hpp"
#include "com/SharedPointer.hpp"
#include "logging/Logger.hpp"

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
  virtual void send(std::string const &itemToSend, int rankReceiver) override;

  /// Sends an array of integer values.
  virtual void send(const int *itemsToSend, int size, int rankReceiver) override;

  /// Asynchronously sends an array of integer values.
  virtual PtrRequest aSend(const int *itemsToSend, int size, int rankReceiver) override;

  /// Sends an array of double values.
  virtual void send(const double *itemsToSend, int size, int rankReceiver) override;

  /// Asynchronously sends an array of double values.
  virtual PtrRequest aSend(const double *itemsToSend, int size, int rankReceiver) override;

  virtual PtrRequest aSend(std::vector<double> const &itemsToSend, int rankReceiver) override;

  /**
   * @brief Sends a double to process with given rank.
   *
   * Default MPI point-to-point communication is used.
   */
  virtual void send(double itemToSend, int rankReceiver) override;

  /// Asynchronously sends a double to process with given rank.
  virtual PtrRequest aSend(const double &itemToSend, int rankReceiver) override;

  /**
   * @brief Sends an int to process with given rank.
   *
   * Default MPI point-to-point communication is used.
   */
  virtual void send(int itemToSend, int rankReceiver) override;

  /// Asynchronously sends an int to process with given rank.
  virtual PtrRequest aSend(const int &itemToSend, int rankReceiver) override;

  /**
   * @brief Sends a bool to process with given rank.
   *
   * Default MPI point-to-point communication is used.
   */
  virtual void send(bool itemToSend, int rankReceiver) override;

  /// Asynchronously sends a bool to process with given rank.
  virtual PtrRequest aSend(const bool &itemToSend, int rankReceiver) override;

  /**
   * @brief Receives a std::string from process with given rank.
   *
   * Default MPI point-to-point communication is used.
   */
  virtual void receive(std::string &itemToReceive, int rankSender) override;

  /// Receives an array of integer values.
  virtual void receive(int *itemsToReceive, int size, int rankSender) override;

  /// Receives an array of double values.
  virtual void receive(double *itemsToReceive, int size, int rankSender) override;

  /// Asynchronously receives an array of double values.
  virtual PtrRequest aReceive(double *itemsToReceive,
                              int     size,
                              int     rankSender) override;

  virtual PtrRequest aReceive(std::vector<double> &itemsToReceive, int rankSender) override;

  /**
   * @brief Receives a double from process with given rank.
   *
   * Default MPI point-to-point communication is used.
   */
  virtual void receive(double &itemToReceive, int rankSender) override;

  /// Asynchronously receives a double from process with given rank.
  virtual PtrRequest aReceive(double &itemToReceive, int rankSender) override;

  /**
   * @brief Receives an int from process with given rank.
   *
   * Default MPI point-to-point communication is used.
   */
  virtual void receive(int &itemToReceive, int rankSender) override;

  /// Asynchronously receives an int from process with given rank.
  virtual PtrRequest aReceive(int &itemToReceive, int rankSender) override;

  /**
   * @brief Receives a bool from process with given rank.
   *
   * Default MPI point-to-point communication is used.
   */
  virtual void receive(bool &itemToReceive, int rankSender) override;

  /// Asynchronously receives a bool from process with given rank.
  virtual PtrRequest aReceive(bool &itemToReceive, int rankSender) override;

  void send(std::vector<int> const &v, int rankReceiver) override;
  void receive(std::vector<int> &v, int rankSender) override;

  void send(std::vector<double> const &v, int rankReceiver) override;
  void receive(std::vector<double> &v, int rankSender) override;

protected:
  /// Returns the communicator.
  virtual MPI_Comm &communicator(int rank) = 0;

  virtual int rank(int rank) = 0;

private:
  logging::Logger _log{"com::MPICommunication"};
};
} // namespace com
} // namespace precice

#endif // not PRECICE_NO_MPI
