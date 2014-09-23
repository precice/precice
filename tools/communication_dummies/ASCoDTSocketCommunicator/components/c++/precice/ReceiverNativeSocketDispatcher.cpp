#include "precice/ReceiverNativeSocketDispatcher.h"
#ifdef _WIN32
#include <winsock2.h>
#include <ws2tcpip.h>
#define bzero(b,len) (memset((b), '\0', (len)), (void) 0)  
#define bcopy(b1,b2,len) (memmove((b2), (b1), (len)), (void) 0)
#else
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <unistd.h>
#include <arpa/inet.h>
#include <fcntl.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <algorithm>
#include <iostream>
#include <string.h>
#include <sstream>

#ifdef _WIN32
#include <windows.h>
DWORD WINAPI Receiver_dispatch_loop(void *arg){
     ((precice::ReceiverNativeSocketDispatcher*) arg)->dispatch();
}
#else

#include <pthread.h>
void* Receiver_dispatch_loop(void *arg){
     ((precice::ReceiverNativeSocketDispatcher*) arg)->dispatch();
}
#endif

precice::ReceiverNativeSocketDispatcher::ReceiverNativeSocketDispatcher(char* host,int port,int bufferSize):
     _buffer_size(bufferSize){
     _rcvBuffer=new char[_buffer_size];
     _sendBuffer=new char[_buffer_size];
     
     precice::ReceiverNativeSocketDispatcher::open(host,port,_sockfd,_newsockfd);
     #ifdef _WIN32
     CreateThread(NULL, 0, Receiver_dispatch_loop, this, 0, NULL);
     #else
     pthread_t task;
     pthread_create(&task, NULL,Receiver_dispatch_loop, this);
	 #endif
}

precice::ReceiverNativeSocketDispatcher::ReceiverNativeSocketDispatcher(int port,int bufferSize):
      _buffer_size(bufferSize){
     _rcvBuffer=new char[_buffer_size];
     _sendBuffer=new char[_buffer_size];
     //precice::ReceiverNativeSocketDispatcher::open(port,_sockfd,_newsockfd);
      #ifdef _WIN32
     CreateThread(NULL, 0, Receiver_dispatch_loop, this, 0, NULL);
     #else
     pthread_t task;
     pthread_create(&task, NULL,Receiver_dispatch_loop, this);
	 #endif
}
precice::ReceiverNativeSocketDispatcher::~ReceiverNativeSocketDispatcher(){
     delete [] _rcvBuffer;
     delete [] _sendBuffer;
     precice::ReceiverNativeSocketDispatcher::close(_sockfd,_newsockfd);
     
}


void precice::ReceiverNativeSocketDispatcher::open(char* hostname,int port,
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
&newsockfd){
          std::stringstream port_str;
		  port_str<<port;
          struct addrinfo *result = NULL,
     	                   hints;
          bzero( &hints, sizeof(hints) );
          hints.ai_family = AF_UNSPEC;
          hints.ai_socktype = SOCK_STREAM;
          hints.ai_protocol = IPPROTO_TCP;
          getaddrinfo(hostname, port_str.str().c_str(), &hints, &result);
          sockfd = socket(result->ai_family, result->ai_socktype, 
            result->ai_protocol);
          assert(sockfd >= 0);
          
          int tries=60;
                while(tries>0){
                 if (::connect(sockfd, result->ai_addr, (int)result->ai_addrlen)==0)
                    break;
                         tries--;
#ifdef _WIN32
						 Sleep(1000);
#else
						 sleep(1);
#endif	
                }   
         newsockfd=sockfd;
        
}

/*
void precice::ReceiverNativeSocketDispatcher::open(int port,int &sockfd,int &newsockfd){
          socklen_t clilen;

          struct sockaddr_in serv_addr, cli_addr;

          sockfd = socket(AF_INET, SOCK_STREAM, 0);
          assert (sockfd>=0);

          bzero((char *) &serv_addr, sizeof(serv_addr));
          serv_addr.sin_family = AF_INET;
          serv_addr.sin_addr.s_addr = INADDR_ANY;
          serv_addr.sin_port = htons(port);
          assert (bind(sockfd, (struct sockaddr *) &serv_addr,
                    sizeof(serv_addr)) >= 0);
          listen(sockfd,5);
          clilen = sizeof(cli_addr);
          newsockfd = accept(sockfd,
                    (struct sockaddr *) &cli_addr,
                    &clilen);
          assert (newsockfd >= 0);
}
*/

void precice::ReceiverNativeSocketDispatcher::sendData(char* data, size_t numberOfBytes, char* sendBuffer,
#ifdef _WIN32
    SOCKET
#else
    int
#endif
newsockfd,int bufferSize){
     char* data_ptr=(char*)data;
     int remaining_bytes_to_send=0,total_send_bytes=0,send_bytes=0,n=0;

     //clear buffer
     bzero(sendBuffer,bufferSize);
     while(total_send_bytes<numberOfBytes){
          remaining_bytes_to_send=(numberOfBytes-total_send_bytes<=bufferSize)?numberOfBytes-total_send_bytes:bufferSize;
          memcpy(sendBuffer,data_ptr,remaining_bytes_to_send);
          send_bytes=0;
          char* send_buffer_ptr=sendBuffer;
          while(send_bytes<bufferSize){
               n = 
#ifdef _WIN32
				   send(
#else

				   write(
#endif
				   newsockfd,send_buffer_ptr,bufferSize-send_bytes
#ifdef _WIN32
				   ,0
#else

#endif
				   );
               if(n>0){
                    send_bytes+=n;
                    send_buffer_ptr+=n;
               }
          }
          total_send_bytes+=send_bytes;
          data_ptr+=send_bytes;
     }

}


void precice::ReceiverNativeSocketDispatcher::readData(char* data,size_t size_of_data,char* readBuffer,
#ifdef _WIN32
    SOCKET
#else
    int
#endif
newsockfd, int bufferSize){
 		  bzero(readBuffer,bufferSize);
          int bytes_received=0;
          int total_bytes_received=0;
          int local_bytes_received=0;
          int bytes_to_copy=0;
          char* data_ptr=data;
          char* buffer_ptr;

          while(total_bytes_received<size_of_data){
               bytes_received=0;
               buffer_ptr=readBuffer;
               while(bytes_received<bufferSize){
                    local_bytes_received = 
#ifdef _WIN32
						recv(
#else
						read(
#endif
							newsockfd,buffer_ptr,bufferSize-bytes_received
#ifdef _WIN32
							,0
#else
#endif
							);
                    if(local_bytes_received>0){
                         bytes_received+=local_bytes_received;
                         buffer_ptr+=local_bytes_received;
                    }
               }
               bytes_to_copy=(total_bytes_received+bytes_received>size_of_data)?size_of_data-total_bytes_received:bytes_received;
               memcpy(data_ptr,readBuffer,bytes_to_copy);
               data_ptr+=bytes_to_copy;
               total_bytes_received+=bytes_to_copy;
          }


}


void precice::ReceiverNativeSocketDispatcher::close(
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
&newsockfd){
#ifdef _WIN32
	 if(newsockfd>=0)
         closesocket(newsockfd);
     if(sockfd>=0)
         closesocket(sockfd);
#else

     if(newsockfd>=0)
         ::close(newsockfd);
     if(sockfd>=0)
         ::close(sockfd);
#endif
}

void connect_Receiver(void* ref,precice::Receiver* port){
    
     ((precice::ReceiverNativeSocketDispatcher*)ref)->connect(port);
}

void disconnect_Receiver(void* ref,precice::Receiver* port){
     ((precice::ReceiverNativeSocketDispatcher*)ref)->disconnect(port);
     delete port;
}

void precice::ReceiverNativeSocketDispatcher::dispatch(){
     int methodId=0;
     long long ref;
     void (*invokers[2])(void*,precice::Receiver*);
     invokers[0]=&connect_Receiver;
     invokers[1]=&disconnect_Receiver;
     do{
          readData((char*)&methodId,sizeof(int),_rcvBuffer,_newsockfd,_buffer_size);
          if(methodId!=-1){
               readData((char*)&ref,sizeof(long long),_rcvBuffer,_newsockfd,_buffer_size);
               (*invokers[methodId])(this,(precice::Receiver*)ref);
          }
     }while(methodId!=-1);
}

