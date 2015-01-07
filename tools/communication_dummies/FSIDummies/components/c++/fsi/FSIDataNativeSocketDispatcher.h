#ifndef FSI_FSIDATANATIVESOCKETDISPATCHER_H_
#define FSI_FSIDATANATIVESOCKETDISPATCHER_H_ 

#include "fsi/FSIDataNativeDispatcher.h"
#include <iostream>
#include <vector>
#ifdef _WIN32
#include <winsock2.h>
#endif
namespace fsi { 

     class FSIDataNativeSocketDispatcher;
}



#ifdef __cplusplus
  extern "C" {
#endif


          


#ifdef __cplusplus
  }
#endif


class fsi::FSIDataNativeSocketDispatcher: public fsi::FSIDataNativeDispatcher{
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
    FSIDataNativeSocketDispatcher(char* host,int port,int bufferSize);
    FSIDataNativeSocketDispatcher(int port,int bufferSize);
    ~FSIDataNativeSocketDispatcher();
    void dispatch();
  
    
   
    
};

#endif
