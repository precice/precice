#ifndef FSI_FSIDATACXX2SOCKETPLAINPORT_H_
#define FSI_FSIDATACXX2SOCKETPLAINPORT_H_ 

#include "fsi/FSIData.h"
#include <iostream>
#include <string>
#ifdef _WIN32
#include <winsock2.h>
#endif
namespace fsi { 

     class FSIDataCxx2SocketPlainPort;
}

class fsi::FSIDataCxx2SocketPlainPort: public fsi::FSIData{
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
    FSIDataCxx2SocketPlainPort(char* host,int port,int buffer_size);
     FSIDataCxx2SocketPlainPort(int port,int buffer_size);
    ~FSIDataCxx2SocketPlainPort();
    //int getSockfd();
    //int getNewsockfd();
    
    void transferData(const int* coordId, const int coordId_len,const double* data, const int data_len);  
    void transferDataParallel(const int* coordId, const int coordId_len,const double* data, const int data_len);
   
    void dataAck(int& ack);  
    void dataAckParallel(int& ack);
   
};

#endif
