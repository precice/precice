
#include <iostream>
#include <string.h>
#include "precice/ReceiverC2CxxSocketPlainPort.h"
#include "precice/ReceiverCxx2SocketPlainPort.h"
#include "precice/Receiver.h"
extern "C" {

#ifdef _WIN32
void PRECICE_RECEIVERC2SOCKET_PLAIN_PORT_CREATE_CLIENT_INSTANCE(long long* ptr,char* host,int& port,int& buffer_size){
#else
void precice_receiverc2socket_plain_port_create_client_instance_(long long* ptr,char* host,int& port,int& buffer_size){
#endif
   *ptr=(long long)new precice::ReceiverCxx2SocketPlainPort(
        host,
        port,
        buffer_size
   );
     

}

#ifdef _WIN32
void PRECICE_RECEIVERC2SOCKET_PLAIN_PORT_CREATE_SERVER_INSTANCE(long long* ptr,int& port,int& buffer_size){
#else
void precice_receiverc2socket_plain_port_create_server_instance_(long long* ptr,int& port,int& buffer_size){
#endif
   *ptr=(long long)new precice::ReceiverCxx2SocketPlainPort(
        port,
        buffer_size
   );
     

}

#ifdef _WIN32
void PRECICE_RECEIVERC2SOCKET_PLAIN_PORT_DESTROY_INSTANCE(long long *ptr){
#else
void precice_receiverc2socket_plain_port_destroy_instance_(long long *ptr){
#endif
     precice::ReceiverCxx2SocketPlainPort* c_ptr = (precice::ReceiverCxx2SocketPlainPort*)*ptr;
     if(c_ptr!=NULL){
         delete c_ptr;
         c_ptr = NULL;
     }
}




#ifdef _WIN32
void PRECICE_RECEIVERC2SOCKET_PLAIN_PORT_RECEIVE(long long* ref,double* data,int* index,int* rank,int* tag){
     
     ((precice::ReceiverCxx2SocketPlainPort*)*ref)->receive(*data,*index,*rank,*tag);
}
#else
void precice_receiverc2socket_plain_port_receive_(long long* ref,double* data,int* index,int* rank,int* tag){
     
     ((precice::ReceiverCxx2SocketPlainPort*)*ref)->receive(*data,*index,*rank,*tag);
}
#endif
}
