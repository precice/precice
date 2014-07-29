#ifndef FSI_FSITESTCXX2SOCKETPLAINPORT_H_
#define FSI_FSITESTCXX2SOCKETPLAINPORT_H_ 

#include "fsi/FSITest.h"
#include <iostream>
#include <string>
#ifdef _WIN32
#include <winsock2.h>
#endif
namespace fsi { 

     class FSITestCxx2SocketPlainPort;
}

class fsi::FSITestCxx2SocketPlainPort: public fsi::FSITest{
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
    FSITestCxx2SocketPlainPort(char* host,int port,int buffer_size);
     FSITestCxx2SocketPlainPort(int port,int buffer_size);
    ~FSITestCxx2SocketPlainPort();
    //int getSockfd();
    //int getNewsockfd();
    
    void test();  
    void testParallel();
   
};

#endif
