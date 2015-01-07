#ifndef PRECICE_MAINNATIVESOCKETDISPATCHER_H_
#define PRECICE_MAINNATIVESOCKETDISPATCHER_H_ 

#include "precice/MainNativeDispatcher.h"
#include <iostream>
#include <vector>
#ifdef _WIN32
#include <winsock2.h>
#endif
namespace precice { 

     class MainNativeSocketDispatcher;
}



#ifdef __cplusplus
  extern "C" {
#endif


          


#ifdef __cplusplus
  }
#endif


class precice::MainNativeSocketDispatcher: public precice::MainNativeDispatcher{
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
    
    void open(char* hostname,int port,
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
    //void open(int port,int &sockfd,int &newsockfd);
    void sendData(char* data, size_t numberOfBytes, char* sendBuffer,
    #ifdef _WIN32
    SOCKET
	#else
    int
	#endif 
    newsockfd,
    int bufferSize);
    void readData(char* data,size_t size_of_data,char* readBuffer,
    #ifdef _WIN32
    SOCKET
	#else
    int
	#endif 
    newsockfd,
    int bufferSize);
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
    MainNativeSocketDispatcher(char* host,int port,int bufferSize);
    MainNativeSocketDispatcher(int port,int bufferSize);
    ~MainNativeSocketDispatcher();
    void dispatch();
  
    
   
    
};

#endif
