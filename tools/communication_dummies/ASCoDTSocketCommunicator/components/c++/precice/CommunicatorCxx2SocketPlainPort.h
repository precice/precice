#ifndef PRECICE_COMMUNICATORCXX2SOCKETPLAINPORT_H_
#define PRECICE_COMMUNICATORCXX2SOCKETPLAINPORT_H_ 

#include "precice/Communicator.h"
#include <iostream>
#include <string>
#ifdef _WIN32
#include <winsock2.h>
#endif
namespace precice { 

     class CommunicatorCxx2SocketPlainPort;
}

class precice::CommunicatorCxx2SocketPlainPort: public precice::Communicator{
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
    void open_client(char* hostname,int port,
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
    CommunicatorCxx2SocketPlainPort(char* host,int port,int buffer_size);
     CommunicatorCxx2SocketPlainPort(int port,int buffer_size);
    ~CommunicatorCxx2SocketPlainPort();
    //int getSockfd();
    //int getNewsockfd();
    
    void setData(const double* data, const int data_len);  
    void setDataParallel(const double* data, const int data_len);
   
};

#endif
