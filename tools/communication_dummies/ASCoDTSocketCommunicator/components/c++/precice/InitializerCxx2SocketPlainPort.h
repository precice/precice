#ifndef PRECICE_INITIALIZERCXX2SOCKETPLAINPORT_H_
#define PRECICE_INITIALIZERCXX2SOCKETPLAINPORT_H_ 

#include "precice/Initializer.h"
#include <iostream>
#include <string>
#ifdef _WIN32
#include <winsock2.h>
#endif
namespace precice { 

     class InitializerCxx2SocketPlainPort;
}

class precice::InitializerCxx2SocketPlainPort: public precice::Initializer{
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
    InitializerCxx2SocketPlainPort(char const* host,int port,int buffer_size);
     InitializerCxx2SocketPlainPort(int port,int buffer_size);
    ~InitializerCxx2SocketPlainPort();
    //int getSockfd();
    //int getNewsockfd();
    
    void acknowledge(const int identifier,int& tag);  
    void acknowledgeParallel(const int identifier,int& tag);
   
    void initialize(const std::string* addresses, const int addresses_len,const int* vertexes, const int vertexes_len);  
    void initializeParallel(const std::string* addresses, const int addresses_len,const int* vertexes, const int vertexes_len);
   
};

#endif
