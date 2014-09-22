
#include <iostream>
#include <string.h>
#include "precice/CommunicatorC2CxxSocketPlainPort.h"
#include "precice/CommunicatorCxx2SocketPlainPort.h"
#include "precice/Communicator.h"
extern "C" {

#ifdef _WIN32
void PRECICE_COMMUNICATORC2SOCKET_PLAIN_PORT_CREATE_CLIENT_INSTANCE(long long* ptr,char* host,int& port,int& buffer_size){
#else
void precice_communicatorc2socket_plain_port_create_client_instance_(long long* ptr,char* host,int& port,int& buffer_size){
#endif
   *ptr=(long long)new precice::CommunicatorCxx2SocketPlainPort(
        host,
        port,
        buffer_size
   );
     

}

#ifdef _WIN32
void PRECICE_COMMUNICATORC2SOCKET_PLAIN_PORT_CREATE_SERVER_INSTANCE(long long* ptr,int& port,int& buffer_size){
#else
void precice_communicatorc2socket_plain_port_create_server_instance_(long long* ptr,int& port,int& buffer_size){
#endif
   *ptr=(long long)new precice::CommunicatorCxx2SocketPlainPort(
        port,
        buffer_size
   );
     

}

#ifdef _WIN32
void PRECICE_COMMUNICATORC2SOCKET_PLAIN_PORT_DESTROY_INSTANCE(long long *ptr){
#else
void precice_communicatorc2socket_plain_port_destroy_instance_(long long *ptr){
#endif
     precice::CommunicatorCxx2SocketPlainPort* c_ptr = (precice::CommunicatorCxx2SocketPlainPort*)*ptr;
     if(c_ptr!=NULL){
         delete c_ptr;
         c_ptr = NULL;
     }
}




#ifdef _WIN32
void PRECICE_COMMUNICATORC2SOCKET_PLAIN_PORT_SETDATA(long long* ref,double* data,int* index,int* rank,int* tag){
     
     ((precice::CommunicatorCxx2SocketPlainPort*)*ref)->setData(*data,*index,*rank,*tag);
}
#else
void precice_communicatorc2socket_plain_port_setdata_(long long* ref,double* data,int* index,int* rank,int* tag){
     
     ((precice::CommunicatorCxx2SocketPlainPort*)*ref)->setData(*data,*index,*rank,*tag);
}
#endif
}
