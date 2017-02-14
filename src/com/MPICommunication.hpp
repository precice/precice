#ifndef PRECICE_NO_MPI

#pragma once

#include <mpi.h>
#include "logging/Logger.hpp"
#include "com/Communication.hpp"

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
  /**
   * @brief Constructor, takes communicator for default communication.
   */
  MPICommunication();

  /**
   * @brief Destructor, empty.
   */
  virtual ~MPICommunication() {
  }

  virtual void
  startSendPackage(int rankReceiver) {
  }

  virtual void
  finishSendPackage() {
  }

  /**
   * @return rankSender
   */
  virtual int
  startReceivePackage(int rankSender) {
    return rankSender;
  }

  virtual void
  finishReceivePackage() {
  }

  /**
   * @brief Sends a std::string to process with given rank.
   *
   * Default MPI point-to-point communication is used.
   */
  virtual void send(std::string const& itemToSend, int rankReceiver);

  /**
   * @brief Sends an array of integer values.
   */
  virtual void send(int* itemsToSend, int size, int rankReceiver);

  /**
   * @brief Asynchronously sends an array of integer values.
   */
  virtual PtrRequest aSend(int* itemsToSend,
                                       int size,
                                       int rankReceiver);

  /**
   * @brief Sends an array of double values.
   */
  virtual void send(double* itemsToSend, int size, int rankReceiver);

  /**
   * @brief Asynchronously sends an array of double values.
   */
  virtual PtrRequest aSend(double* itemsToSend,
                                       int size,
                                       int rankReceiver);

  /**
   * @brief Sends a double to process with given rank.
   *
   * Default MPI point-to-point communication is used.
   */
  virtual void send(double itemToSend, int rankReceiver);

  /**
   * @brief Asynchronously sends a double to process with given rank.
   */
  virtual PtrRequest aSend(double* itemToSend, int rankReceiver);

  /**
   * @brief Sends an int to process with given rank.
   *
   * Default MPI point-to-point communication is used.
   */
  virtual void send(int itemToSend, int rankReceiver);

  /**
   * @brief Asynchronously sends an int to process with given rank.
   */
  virtual PtrRequest aSend(int* itemToSend, int rankReceiver);

  /**
   * @brief Sends a bool to process with given rank.
   *
   * Default MPI point-to-point communication is used.
   */
  virtual void send(bool itemToSend, int rankReceiver);

  /**
   * @brief Asynchronously sends a bool to process with given rank.
   */
  virtual PtrRequest aSend(bool* itemToSend, int rankReceiver);

  /**
   * @brief Receives a std::string from process with given rank.
   *
   * Default MPI point-to-point communication is used.
   */
  virtual void receive(std::string& itemToReceive, int rankSender);

  /**
   * @brief Receives an array of integer values.
   */
  virtual void receive(int* itemsToReceive, int size, int rankSender);

  /**
   * @brief Asynchronously receives an array of integer values.
   */
  virtual PtrRequest aReceive(int* itemsToReceive,
                                          int size,
                                          int rankSender);

  /**
   * @brief Receives an array of double values.
   */
  virtual void receive(double* itemsToReceive, int size, int rankSender);

  /**
   * @brief Asynchronously receives an array of double values.
   */
  virtual PtrRequest aReceive(double* itemsToReceive,
                                          int size,
                                          int rankSender);

  /**
   * @brief Receives a double from process with given rank.
   *
   * Default MPI point-to-point communication is used.
   */
  virtual void receive(double& itemToReceive, int rankSender);

  /**
   * @brief Asynchronously receives a double from process with given rank.
   */
  virtual PtrRequest aReceive(double* itemToReceive,
                                          int rankSender);

  /**
   * @brief Receives an int from process with given rank.
   *
   * Default MPI point-to-point communication is used.
   */
  virtual void receive(int& itemToReceive, int rankSender);

  /**
   * @brief Asynchronously receives an int from process with given rank.
   */
  virtual PtrRequest aReceive(int* itemToReceive, int rankSender);

  /**
   * @brief Receives a bool from process with given rank.
   *
   * Default MPI point-to-point communication is used.
   */
  virtual void receive(bool& itemToReceive, int rankSender);

  /// Asynchronously receives a bool from process with given rank.
  virtual PtrRequest aReceive(bool* itemToReceive, int rankSender);

protected:
  /// Returns the communicator.
  virtual MPI_Comm& communicator(int rank) = 0;

  virtual int rank(int rank) = 0;

private:

  static logging::Logger _log;
};

}} // namespace precice, com


#endif // not PRECICE_NO_MPI
