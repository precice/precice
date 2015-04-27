// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at
// http://www5.in.tum.de/wiki/index.php/PreCICE_License

#ifndef PRECICE_COM_FILE_COMMUNICATION_HPP_
#define PRECICE_COM_FILE_COMMUNICATION_HPP_

#include "Communication.hpp"

#include "tarch/logging/Log.h"

#include <fstream>
#include <map>
#include <string>

namespace precice {
namespace com {
/**
 * @brief Implementation of Communication interface by using files.
 *
 * NOTE:
 * Asynchronous sending methods are not implemented.
 */
class FileCommunication : public Communication {
public:
  /**
   * @brief Constructor.
   */
  FileCommunication(bool binaryMode, const std::string& communicationDirectory);

  /**
   * @brief Destructor, empty.
   */
  virtual ~FileCommunication();

  /**
   * @brief Returns true, if a connection to a remote participant has been
   * setup.
   */
  virtual bool isConnected();

  /**
   * @brief Returns 1, since file com. is not prepared for parallelity.
   *
   * Precondition: a connection to the remote participant has been setup.
   */
  virtual int
  getRemoteCommunicatorSize() {
    return 1;
  };

  /**
   * @brief Gathers information about files to write and read to/from.
   */
  virtual void acceptConnection(std::string const& nameAcceptor,
                                std::string const& nameRequester,
                                int acceptorProcessRank,
                                int acceptorCommunicatorSize);

  virtual void acceptConnectionAsServer(std::string const& nameAcceptor,
                                        std::string const& nameRequester,
                                        int requesterCommunicatorSize);

  /**
   * @brief Gathers information about files to write and read to/from.
   */
  virtual void requestConnection(std::string const& nameAcceptor,
                                 std::string const& nameRequester,
                                 int requesterProcessRank,
                                 int requesterCommunicatorSize);

  virtual int
  requestConnectionAsClient(std::string const& nameAcceptor,
                            std::string const& nameRequester);

  /**
   * @brief Doesn't do anything here.
   */
  virtual void closeConnection();

  /**
   * @brief Opens a file to write all send messages into it.
   */
  virtual void startSendPackage(int rankReceiver);

  /**
   * @brief Closes the file with written messages and makes it available.
   */
  virtual void finishSendPackage();

  /**
   * @brief Opens the file that holds all messages to receive.
   */
  virtual int startReceivePackage(int rankSender);

  /**
   * @brief Closes and removes the file with read messages.
   */
  virtual void finishReceivePackage();

  /**
   * @brief Sends a std::string to process with given rank.
   */
  virtual void send(std::string const& itemToSend, int rankReceiver);

  /**
   * @brief Sends an array of integer values.
   */
  virtual void send(int* itemsToSend, int size, int rankReceiver);

  /**
   * @brief Asynchronously sends an array of integer values.
   */
  virtual Request::SharedPointer
  aSend(int* itemsToSend, int size, int rankReceiver);

  /**
   * @brief Sends an array of double values.
   */
  virtual void send(double* itemsToSend, int size, int rankReceiver);

  /**
   * @brief Asynchronously sends an array of double values.
   */
  virtual Request::SharedPointer
  aSend(double* itemsToSend, int size, int rankReceiver);

  /**
   * @brief Sends a double to process with given rank.
   */
  virtual void send(double itemToSend, int rankReceiver);

  /**
   * @brief Asynchronously sends a double to process with given rank.
   */
  virtual Request::SharedPointer
  aSend(double* itemToSend, int rankReceiver);

  /**
   * @brief Sends an int to process with given rank.
   */
  virtual void send(int itemToSend, int rankReceiver);

  /**
   * @brief Asynchronously sends an int to process with given rank.
   */
  virtual Request::SharedPointer
  aSend(int* itemToSend, int rankReceiver);

  /**
   * @brief Sends a bool to process with given rank.
   */
  virtual void send(bool itemToSend, int rankReceiver);

  /**
   * @brief Asynchronously sends a bool to process with given rank.
   */
  virtual Request::SharedPointer
  aSend(bool* itemToSend, int rankReceiver);

  /**
   * @brief Receives a std::string from process with given rank.
   */
  virtual int receive(std::string& itemToReceive, int rankSender);

  /**
   * @brief Receives an array of integer values.
   */
  virtual int receive(int* itemsToReceive, int size, int rankSender);

  /**
   * @brief Receives an array of double values.
   */
  virtual int receive(double* itemsToReceive, int size, int rankSender);

  /**
   * @brief Receives a double from process with given rank.
   */
  virtual int receive(double& itemToReceive, int rankSender);

  /**
   * @brief Receives an int from process with given rank.
   */
  virtual int receive(int& itemToReceive, int rankSender);

  /**
   * @brief Receives a bool from process with given rank.
   */
  virtual int receive(bool& itemToReceive, int rankSender);

private:
  // @brief Logging device.
  static tarch::logging::Log _log;

  const int TYPE_DOUBLE;
  const int TYPE_INT;
  const int TYPE_STRING;
  const int TYPE_DOUBLE_VECTOR;
  const int TYPE_INT_VECTOR;
  const int TYPE_BOOL;

  // @brief Map from receiver indices to count of sent messages to that
  // receiver.
  // Used to name files.
  std::map<int, int> _sendIndices;

  // @brief Map from sender indices to count of received messages from that
  // sender.
  // Used to name files.
  std::map<int, int> _receiveIndices;

  int _currentPackageRank;

  // @brief Set on communication setup, and used when receiving from any sender.
  // int _remoteCommunicatorSize;

  int _localRank;

  std::string _nameLocal;

  std::string _nameRemote;

  std::ofstream _sendFile;

  std::ifstream _receiveFile;

  bool _isConnected;

  std::ios_base::openmode _sendmode;

  std::ios_base::openmode _receivemode;

  bool _binary;

  // @brief Files for communication are written to that directory.
  std::string _comDirectory;

  /**
   * @brief Renames a file with send data, such that it is available to be read.
   */
  void makeSendFileAvailable(int rank, int index);

  /**
   * @brief Renames a file with read data, such that it can be read privately.
   */
  void makeReceiveFileUnavailable(int rank, int index);

  /**
   * @brief Removes a read (hidden) file.
   */
  void removeReceiveFile(int rank, int index);

  /**
   * @brief Returns a filename for send data for the given parameters.
   *
   * @param hidden [IN] A hidden send file is not ready for reception.
   * @param rankRemote [IN] Rank number of remote process reading the file.
   * @param index [IN] Counter of send files, must match on both sides.
   */
  std::string getSendFilename(bool hidden, int rankRemote, int index);

  /**
   * @brief Returns a filename for receive data for the given parameters.
   *
   * @param hidden [IN] A hidden receive file is (ready for) being read.
   * @param rankRemote [IN] Rank number of remote process writing the file.
   * @param index [IN] Counter of send files, must match on both sides.
   */
  std::string getReceiveFilename(bool hidden, int rankRemote, int index);
};
}
} // namespace precice, com

#endif /* PRECICE_COM_FILE_COMMUNICATION_HPP_ */
