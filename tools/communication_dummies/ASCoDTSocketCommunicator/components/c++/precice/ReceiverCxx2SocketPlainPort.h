#ifndef PRECICE_RECEIVERCXX2SOCKETPLAINPORT_H_
#define PRECICE_RECEIVERCXX2SOCKETPLAINPORT_H_ 

#include "precice/Receiver.h"
#include <iostream>
#include <string>
#ifdef _WIN32
#include <winsock2.h>
#endif
namespace precice { 

     class ReceiverCxx2SocketPlainPort;
}

class precice::ReceiverCxx2SocketPlainPort: public precice::Receiver{
  private:
    #ifdef _WIN32
    SOCKET
    #else
    int
    #endif 
    _sockfd;
    #ifdef _WIN32
    SOCKET
    #else
    int
    #endif
    _newsockfd;
    int _buffer_size;
    char *_rcvBuffer;
    char *_sendBuffer;
    void open_client(char const* hostname,int port,
    #ifdef _WIN32
    SOCKET
    #else
    int
    #endif 
    &sockfd,
    #ifdef _WIN32
    SOCKET
    #else
    int
    #endif 
    &newsockfd);
    //void open_server(int port,int &sockfd,int &newsockfd);
    void sendData(char* data, size_t numberOfBytes, char* sendBuffer,
    #ifdef _WIN32
    SOCKET
    #else
    int
    #endif 
    newsockfd,int bufferSize);
    void readData(char* data,size_t size_of_data,char* readBuffer,
    #ifdef _WIN32
    SOCKET
    #else
    int
    #endif 
    newsockfd, int bufferSize);
    void close(
    #ifdef _WIN32
    SOCKET
    #else
    int
    #endif 
    &sockfd,
    #ifdef _WIN32
    SOCKET
    #else
    int
    #endif 
    &newsockfd);
  public:
    ReceiverCxx2SocketPlainPort(char const* host,int port,int buffer_size);
     ReceiverCxx2SocketPlainPort(int port,int buffer_size);
    ~ReceiverCxx2SocketPlainPort();
    //int getSockfd();
    //int getNewsockfd();
    
    void receive(const double data,const int index,const int rank,int& tag);  
    void receiveParallel(const double data,const int index,const int rank,int& tag);
   
};

#endif
