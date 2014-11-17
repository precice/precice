
#include <iostream>
#include <string.h>
#include "precice/InitializerC2CxxSocketPlainPort.h"
#include "precice/InitializerCxx2SocketPlainPort.h"
#include "precice/Initializer.h"
extern "C" {

#ifdef _WIN32
void PRECICE_INITIALIZERC2SOCKET_PLAIN_PORT_CREATE_CLIENT_INSTANCE(long long* ptr,char* host,int& port,int& buffer_size){
#else
void precice_initializerc2socket_plain_port_create_client_instance_(long long* ptr,char* host,int& port,int& buffer_size){
#endif
   *ptr=(long long)new precice::InitializerCxx2SocketPlainPort(
        host,
        port,
        buffer_size
   );
     

}

#ifdef _WIN32
void PRECICE_INITIALIZERC2SOCKET_PLAIN_PORT_CREATE_SERVER_INSTANCE(long long* ptr,int& port,int& buffer_size){
#else
void precice_initializerc2socket_plain_port_create_server_instance_(long long* ptr,int& port,int& buffer_size){
#endif
   *ptr=(long long)new precice::InitializerCxx2SocketPlainPort(
        port,
        buffer_size
   );
     

}

#ifdef _WIN32
void PRECICE_INITIALIZERC2SOCKET_PLAIN_PORT_DESTROY_INSTANCE(long long *ptr){
#else
void precice_initializerc2socket_plain_port_destroy_instance_(long long *ptr){
#endif
     precice::InitializerCxx2SocketPlainPort* c_ptr = (precice::InitializerCxx2SocketPlainPort*)*ptr;
     if(c_ptr!=NULL){
         delete c_ptr;
         c_ptr = NULL;
     }
}




#ifdef _WIN32
void PRECICE_INITIALIZERC2SOCKET_PLAIN_PORT_ACKNOWLEDGE(long long* ref,int* identifier,int* tag){
     
     ((precice::InitializerCxx2SocketPlainPort*)*ref)->acknowledge(*identifier,*tag);
}
#else
void precice_initializerc2socket_plain_port_acknowledge_(long long* ref,int* identifier,int* tag){
     
     ((precice::InitializerCxx2SocketPlainPort*)*ref)->acknowledge(*identifier,*tag);
}
#endif
#ifdef _WIN32
void PRECICE_INITIALIZERC2SOCKET_PLAIN_PORT_INITIALIZE(long long* ref,char** addresses,int* addresses_len,int* vertexes,int* vertexes_len){
     std::string* addresses_str=new std::string[*addresses_len];
for(int i=0;i<*addresses_len;i++)
addresses_str[i]=addresses[i];

     ((precice::InitializerCxx2SocketPlainPort*)*ref)->initialize(addresses_str,*addresses_len,vertexes,*vertexes_len);
}
#else
void precice_initializerc2socket_plain_port_initialize_(long long* ref,char** addresses,int* addresses_len,int* vertexes,int* vertexes_len){
     std::string* addresses_str=new std::string[*addresses_len];
for(int i=0;i<*addresses_len;i++)
addresses_str[i]=addresses[i];

     ((precice::InitializerCxx2SocketPlainPort*)*ref)->initialize(addresses_str,*addresses_len,vertexes,*vertexes_len);
}
#endif
}
