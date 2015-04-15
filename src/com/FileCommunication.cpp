// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License

#include "FileCommunication.hpp"

#include "utils/Globals.hpp"
#include "utils/Helpers.hpp"

#include <iomanip>
#include <stdio.h>

namespace precice {
namespace com {

tarch::logging::Log FileCommunication:: _log ( "precice::com::FileCommunication" );

FileCommunication:: FileCommunication
(
  bool               binaryMode,
  const std::string& communicationDirectory)
:
  TYPE_DOUBLE ( 0 ),
  TYPE_INT ( 1 ),
  TYPE_STRING ( 2 ),
  TYPE_DOUBLE_VECTOR ( 3 ),
  TYPE_INT_VECTOR ( 4 ),
  TYPE_BOOL ( 5 ),
  _sendIndices (),
  _receiveIndices (),
  _currentPackageRank ( -1 ),
  _localRank (-1),
  _nameLocal (),
  _nameRemote (),
  _sendFile (),
  _receiveFile (),
  _isConnected ( false ),
  _sendmode ( std::ios::trunc|std::ios::out ),
  _receivemode ( std::ios::in ),
  _binary ( binaryMode ),
  _comDirectory(communicationDirectory)
{
  _sendFile.setf ( std::ios::showpoint );
  _sendFile.setf ( std::ios::fixed );
  _sendFile << std::setprecision(16);
  if ( _binary ){
    _sendmode |= std::ios::binary;
    _receivemode |= std::ios::binary;
  }
}

FileCommunication::~FileCommunication() {
  preciceTrace1("~FileCommunication()", _isConnected);

  closeConnection();
}

bool FileCommunication:: isConnected()
{
  return _isConnected;
}

void FileCommunication:: acceptConnection
(
  const std::string& nameAcceptor,
  const std::string& nameRequester,
  int                acceptorProcessRank,
  int                acceptorCommunicatorSize )
{
  preciceTrace2 ( "acceptConnection()", nameAcceptor, nameRequester );
  _nameLocal = nameAcceptor;
  _nameRemote = nameRequester;
  assertion2 ( (acceptorProcessRank >= 0) && (acceptorProcessRank <= acceptorCommunicatorSize),
               acceptorProcessRank, acceptorCommunicatorSize );
  _localRank = acceptorProcessRank;
  remove ( getReceiveFilename(true,0,0).c_str() ); // To cleanup
  remove ( getSendFilename(false,0,0).c_str() );   // To cleanup
  _isConnected = true;
}

void
FileCommunication::acceptConnectionAsServer(std::string const& nameAcceptor,
                                            std::string const& nameRequester,
                                            int requesterCommunicatorSize) {
  preciceError("acceptConnectionAsServer()", "Not implemented!");
}

void FileCommunication:: requestConnection
(
  const std::string& nameAcceptor,
  const std::string& nameRequester,
  int                requesterProcessRank,
  int                requesterCommunicatorSize )
{
  preciceTrace2 ( "requestConnection()", nameAcceptor, nameRequester );
  _nameLocal = nameRequester;
  _nameRemote = nameAcceptor;
  assertion2 ( (requesterProcessRank >= 0) && (requesterProcessRank <= requesterCommunicatorSize),
                requesterProcessRank, requesterCommunicatorSize );
  _localRank = requesterProcessRank;
  remove ( getSendFilename(false,0,0).c_str() ); // To cleanup
  remove ( getReceiveFilename(true,0,0).c_str() ); // To cleanup
  _isConnected = true;
}

int
FileCommunication::requestConnectionAsClient(std::string const& nameAcceptor,
                                             std::string const& nameRequester) {
  preciceError("requestConnectionAsClient()", "Not implemented!");
}

void FileCommunication:: closeConnection()
{
  preciceTrace ( "closeConnection()" );

  if (not isConnected())
    return;

  _isConnected = false;
}

void FileCommunication:: startSendPackage
(
  int rankReceiver )
{
  preciceTrace ( "startSendPackage()" );
  assertion1 ( _currentPackageRank == -1, _currentPackageRank );
  assertion1 ( rankReceiver >= 0, rankReceiver );
  assertion ( rankReceiver != ANY_SENDER );
  _currentPackageRank = rankReceiver;
  int sendIndex;
  std::map<int,int>::iterator iter  =_sendIndices.find(rankReceiver);
  if ( iter == _sendIndices.end() ){
    sendIndex = 1;
    _sendIndices.insert ( std::make_pair(rankReceiver, 1) );
  }
  else {
    iter->second++;
    sendIndex = iter->second;
  }
  std::string filename = getSendFilename (true, rankReceiver, sendIndex);
  _sendFile.open ( filename.c_str(), _sendmode );
  preciceCheck ( _sendFile.is_open(), "startSendPackage()", "Could not open file "
                 << "\"" << &_sendFile << "\" to start send package!" );
}

void FileCommunication:: finishSendPackage()
{
  preciceTrace ( "finishSendPackage()" );
  assertion ( _sendFile.is_open() );
  assertion ( _currentPackageRank != -1 );
  _sendFile.close ();
  assertion1 ( utils::contained(_currentPackageRank, _sendIndices),
               _currentPackageRank );
  int sendIndex = _sendIndices[_currentPackageRank];
  makeSendFileAvailable ( _currentPackageRank, sendIndex );
  _currentPackageRank = -1;
}

int FileCommunication:: startReceivePackage
(
  int rankSender )
{
  preciceTrace ( "startReceivePackage()" );
  assertion1 ( _currentPackageRank == -1, _currentPackageRank );
  assertion1 ( rankSender >= 0, rankSender );
  assertion ( rankSender != ANY_SENDER );
  _currentPackageRank = rankSender;
  int receiveIndex;
  std::map<int,int>::iterator iter  =_receiveIndices.find(rankSender);
  if ( iter == _receiveIndices.end() ){
    receiveIndex = 1;
    _receiveIndices.insert ( std::make_pair(rankSender, 1) );
  }
  else {
    iter->second++;
    receiveIndex = iter->second;
  }
  makeReceiveFileUnavailable ( _currentPackageRank, receiveIndex );
  std::string filename = getReceiveFilename ( true, rankSender, receiveIndex );
  _receiveFile.open ( filename.c_str(), _receivemode );
  preciceCheck ( _receiveFile.is_open(), "startReceivePackage()", "Could not open file "
                 << "\"" << &_receiveFile << "\" to start receive package!" );
  return rankSender;
}

void FileCommunication:: finishReceivePackage()
{
  preciceTrace ( "finishReceivePackage()" );
  assertion ( _receiveFile.is_open() );
  assertion ( _currentPackageRank != -1 );
  _receiveFile.close ();
  assertion1 ( utils::contained(_currentPackageRank, _receiveIndices),
               _currentPackageRank );
  int receiveIndex = _receiveIndices[_currentPackageRank];
  removeReceiveFile ( _currentPackageRank, receiveIndex );
  _currentPackageRank = -1;
}

void FileCommunication:: send
(
  const std::string& itemToSend,
  int                rankReceiver )
{
  preciceTrace ( "send(string)" );
  assertion ( _sendFile.is_open() );
  assertion ( rankReceiver != ANY_SENDER );
  _sendFile.write ( (char*)&TYPE_STRING, sizeof(int) );
  int size = itemToSend.size() + 1;
  _sendFile.write ( (char*)&size, sizeof(int) );
  _sendFile.write ( itemToSend.c_str(), size );
}

void FileCommunication:: send
(
  int* itemsToSend,
  int  size,
  int  rankReceiver )
{
  preciceTrace ( "send(int*)" );
  assertion ( _sendFile.is_open() );
  assertion ( rankReceiver != ANY_SENDER );
  _sendFile.write ( (char*)&TYPE_INT_VECTOR, sizeof(int) );
  _sendFile.write ( (char*)&size, sizeof(int) );
  _sendFile.write ( (char*)itemsToSend, sizeof(int)*size );
}

Request::SharedPointer
FileCommunication::aSend(int* itemsToSend, int size, int rankReceiver) {
  preciceError("aSend()", "Not implemented!");
}

void FileCommunication:: send (
  double* itemsToSend,
  int     size,
  int     rankReceiver )
{
  preciceTrace ( "send(double*)" );
  assertion ( _sendFile.is_open() );
  assertion ( rankReceiver != ANY_SENDER );
  _sendFile.write ( (char*)&TYPE_DOUBLE_VECTOR, sizeof(int) );
  _sendFile.write ( (char*)&size, sizeof(int) );
  _sendFile.write ( (char*)itemsToSend, sizeof(double)*size );
}

Request::SharedPointer
FileCommunication::aSend(double* itemsToSend, int size, int rankReceiver) {
  preciceError("aSend()", "Not implemented!");
}

void FileCommunication:: send
(
  double itemToSend,
  int    rankReceiver )
{
  assertion ( _sendFile.is_open() );
  assertion ( rankReceiver != ANY_SENDER );
  _sendFile.write ( (char*)&TYPE_DOUBLE, sizeof(int) );
  _sendFile.write ( (char*)&itemToSend, sizeof(double) );
}

Request::SharedPointer
FileCommunication::aSend(double* itemToSend, int rankReceiver) {
  preciceError("aSend()", "Not implemented!");
}

void FileCommunication:: send
(
  int itemToSend,
  int rankReceiver )
{
  assertion ( _sendFile.is_open() );
  assertion ( rankReceiver != ANY_SENDER );
  _sendFile.write ( (char*)&TYPE_INT, sizeof(int) );
  _sendFile.write ( (char*)&itemToSend, sizeof(int) );
}

Request::SharedPointer
FileCommunication::aSend(int* itemToSend, int rankReceiver) {
  preciceError("aSend()", "Not implemented!");
}

void FileCommunication:: send
(
  bool itemToSend,
  int  rankReceiver )
{
  assertion ( _sendFile.is_open() );
  assertion ( rankReceiver != ANY_SENDER );
  _sendFile.write ( (char*)&TYPE_BOOL, sizeof(int) );
  _sendFile.write ( (char*)&itemToSend, 1 );
}

Request::SharedPointer
FileCommunication::aSend(bool* itemToSend, int rankReceiver) {
  preciceError("aSend()", "Not implemented!");
}

int FileCommunication:: receive
(
  std::string& itemToReceive,
  int          rankSender )
{
  preciceTrace ( "receive(string)" );
  assertion ( _receiveFile.is_open() );
  assertion ( rankSender != ANY_SENDER );
  int type;
  _receiveFile.read ( (char*)&type, sizeof(int) );
  preciceCheck ( type == TYPE_STRING, "receive(string)",
                 "Receive type is different than string!" );
  int size = 0;
  _receiveFile.read ( (char*)&size, sizeof(int) );
  preciceDebug ( "Size = " << size );
  assertion1 ( size < 500, size );
  char* message = new char[size];
  _receiveFile.read ( message, size );
  itemToReceive = message;
  delete[] message;
  return rankSender;
}

int FileCommunication:: receive
(
  int* itemsToReceive,
  int  size,
  int  rankSender )
{
  preciceTrace2 ( "receive(int*)", size, rankSender );
  assertion ( _receiveFile.is_open() );
  assertion ( rankSender != ANY_SENDER );
  int type;
  _receiveFile.read ( (char*)&type, sizeof(int) );
  preciceCheck ( type == TYPE_INT_VECTOR, "receive(int*)",
                 "Receive type is different than int*!" );
  int writtenSize = 0;
  _receiveFile.read ( (char*)&writtenSize, sizeof(int) );
  assertion2 ( size == writtenSize, size, writtenSize );
  _receiveFile.read ( (char*)itemsToReceive, sizeof(int)*size );
  return rankSender;
}

int FileCommunication:: receive
(
  double* itemsToReceive,
  int     size,
  int     rankSender )
{
  preciceTrace2 ( "receive(double*)", size, rankSender );
  assertion ( _receiveFile.is_open() );
  assertion ( rankSender != ANY_SENDER );
  int type;
  _receiveFile.read ( (char*)&type, sizeof(int) );
  preciceCheck ( type == TYPE_DOUBLE_VECTOR, "receive(double*)",
                 "Receive type is different than double*!" );
  int writtenSize = 0;
  _receiveFile.read ( (char*)&writtenSize, sizeof(int) );
  assertion2 ( size == writtenSize, size, writtenSize );
  _receiveFile.read ( (char*)itemsToReceive, sizeof(double)*size );
  return rankSender;
}

int FileCommunication:: receive
(
   double& itemToReceive,
   int     rankSender )
{
  preciceTrace1 ( "receive(double)", rankSender );
  assertion ( _receiveFile.is_open() );
  assertion ( rankSender != ANY_SENDER );
  int type;
  _receiveFile.read ( (char*)&type, sizeof(int) );
  preciceCheck ( type == TYPE_DOUBLE, "receive(double)",
                 "Receive type is different than double!" );
  _receiveFile.read ( (char*)&itemToReceive, sizeof(double) );
  return rankSender;
}

int FileCommunication:: receive
(
  int& itemToReceive,
  int  rankSender )
{
  preciceTrace1 ( "receive(int)", rankSender );
  assertion ( _receiveFile.is_open() );
  assertion ( rankSender != ANY_SENDER );
  int type;
  _receiveFile.read ( (char*)&type, sizeof(int) );
  preciceCheck ( type == TYPE_INT, "receive(int)",
                 "Receive type is different than int!" );
  _receiveFile.read ( (char*)&itemToReceive, sizeof(int) );
  return rankSender;
}

int FileCommunication:: receive
(
  bool& itemToReceive,
  int   rankSender )
{
  preciceTrace1 ( "receive(bool)", rankSender );
  assertion ( _receiveFile.is_open() );
  assertion ( rankSender != ANY_SENDER );
  int type;
  _receiveFile.read ( (char*)&type, sizeof(int) );
  preciceCheck ( type == TYPE_BOOL, "receive(bool)",
                 "Receive type is different than bool!" );
  _receiveFile.read ( (char*)&itemToReceive, sizeof(bool) );
  return rankSender;
}

void FileCommunication:: makeSendFileAvailable
(
  int rank,
  int index )
{
  preciceTrace2 ( "makeSendFileAvailable()", rank, index );
  std::string filenameHidden = getSendFilename (true, rank, index);
  std::string filename = getSendFilename (false, rank, index);
  if (rename(filenameHidden.c_str(), filename.c_str()) != 0){
    perror ( "Error renaming file:" );
    preciceError ( "finishSendPackage()", "Could not make send file \""
                   << filenameHidden << "\" available!" );
  }
}

void FileCommunication:: makeReceiveFileUnavailable
(
  int rank,
  int index )
{
  preciceTrace2 ( "makeReceiveFileUnavailable()", rank, index );
  std::string filenameHidden = getReceiveFilename (true, rank, index);
  std::string filename = getReceiveFilename (false, rank, index);
  preciceDebug ( "Rename file \"" << filename << "\" to \"" << filenameHidden
                 << "\"" );
  while ( rename(filename.c_str(), filenameHidden.c_str()) != 0 ){
  }
}

void FileCommunication:: removeReceiveFile
(
  int rank,
  int index )
{
  preciceTrace2 ( "removeReceiveFile()", rank, index );
  if ( remove(getReceiveFilename(true,rank,index).c_str() ) != 0){
    preciceError ( "removeReceiveFile()", "Could not remove receive file \""
                   << getReceiveFilename(true,rank,index) << "\"!" );
  }
}

std::string FileCommunication:: getSendFilename
(
  bool hidden,
  int  rankRemote,
  int  index )
{
  std::ostringstream filename;
  if ( hidden ){
    filename << ".hidden.";
  }
  filename << _comDirectory << "precice_com_file.from_" << _nameLocal << ".rank_" << _localRank
           << ".to_" << _nameRemote << ".rank_" << rankRemote << ".index_" << index;
  return filename.str();
}

std::string FileCommunication:: getReceiveFilename
(
  bool hidden,
  int  rankRemote,
  int  index )
{
  std::ostringstream filename;
  if ( hidden ){
    filename << ".hidden.";
  }
  filename << _comDirectory << "precice_com_file" << ".from_" << _nameRemote << ".rank_" << rankRemote
           << ".to_" << _nameLocal << ".rank_" << _localRank << ".index_" << index;
  return filename.str();
}

}} // namespace precice, com
