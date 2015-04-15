// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at
// http://www5.in.tum.de/wiki/index.php/PreCICE_License

#ifndef PRECICE_NO_MPI

#ifndef PRECICE_COM_MPI_COMMUNICATION_HPP_
#define PRECICE_COM_MPI_COMMUNICATION_HPP_

#include "Communication.hpp"

#include "mesh/Data.hpp"
#include "tarch/logging/Log.h"
#include "utils/Dimensions.hpp"

#include <mpi.h>

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
  virtual Request::SharedPointer aSend(int* itemsToSend,
                                       int size,
                                       int rankReceiver);

  /**
   * @brief Sends an array of double values.
   */
  virtual void send(double* itemsToSend, int size, int rankReceiver);

  /**
   * @brief Asynchronously sends an array of double values.
   */
  virtual Request::SharedPointer aSend(double* itemsToSend,
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
  virtual Request::SharedPointer aSend(double* itemToSend, int rankReceiver);

  /**
   * @brief Sends an int to process with given rank.
   *
   * Default MPI point-to-point communication is used.
   */
  virtual void send(int itemToSend, int rankReceiver);

  /**
   * @brief Asynchronously sends an int to process with given rank.
   */
  virtual Request::SharedPointer aSend(int* itemToSend, int rankReceiver);

  /**
   * @brief Sends a bool to process with given rank.
   *
   * Default MPI point-to-point communication is used.
   */
  virtual void send(bool itemToSend, int rankReceiver);

  /**
   * @brief Asynchronously sends a bool to process with given rank.
   */
  virtual Request::SharedPointer aSend(bool* itemToSend, int rankReceiver);

  /**
   * @brief Receives a std::string from process with given rank.
   *
   * Default MPI point-to-point communication is used.
   *
   * @return Rank of sender, which is useful when ANY_SENDER is used.
   */
  virtual int receive(std::string& itemToReceive, int rankSender);

  /**
   * @brief Receives an array of integer values.
   *
   * @return Rank of sender, which is useful when ANY_SENDER is used.
   */
  virtual int receive(int* itemsToReceive, int size, int rankSender);

  /**
   * @brief Receives an array of double values.
   */
  virtual int receive(double* itemsToReceive, int size, int rankSender);

  /**
   * @brief Receives a double from process with given rank.
   *
   * Default MPI point-to-point communication is used.
   *
   * @return Rank of sender, which is useful when ANY_SENDER is used.
   */
  virtual int receive(double& itemToReceive, int rankSender);

  /**
   * @brief Receives an int from process with given rank.
   *
   * Default MPI point-to-point communication is used.
   *
   * @return Rank of sender, which is useful when ANY_SENDER is used.
   */
  virtual int receive(int& itemToReceive, int rankSender);

  /**
   * @brief Receives a bool from process with given rank.
   *
   * Default MPI point-to-point communication is used.
   *
   * @return Rank of sender, which is useful when ANY_SENDER is used.
   */
  virtual int receive(bool& itemToReceive, int rankSender);

protected:
  /**
   * @brief Returns the communicator.
   */
  virtual MPI_Comm& communicator(int rank) = 0;

  virtual int rank(int rank) = 0;

private:
  // @brief Logging device.
  static tarch::logging::Log _log;
};
}
} // namespace precice, com

#endif /* PRECICE_COM_COMMUNICATIONMPI_HPP_ */

#endif // not PRECICE_NO_MPI
