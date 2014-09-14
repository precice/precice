
#include <iostream>
#include <string.h>
#include "precice/MainC2CxxSocketPlainPort.h"
#include "precice/MainCxx2SocketPlainPort.h"
#include "precice/Main.h"
extern "C" {

#ifdef _WIN32
void PRECICE_MAINC2SOCKET_PLAIN_PORT_CREATE_CLIENT_INSTANCE(long long* ptr,char* host,int& port,int& buffer_size){
#else
void precice_mainc2socket_plain_port_create_client_instance_(long long* ptr,char* host,int& port,int& buffer_size){
#endif
   *ptr=(long long)new precice::MainCxx2SocketPlainPort(
        host,
        port,
        buffer_size
   );
     

}

#ifdef _WIN32
void PRECICE_MAINC2SOCKET_PLAIN_PORT_CREATE_SERVER_INSTANCE(long long* ptr,int& port,int& buffer_size){
#else
void precice_mainc2socket_plain_port_create_server_instance_(long long* ptr,int& port,int& buffer_size){
#endif
   *ptr=(long long)new precice::MainCxx2SocketPlainPort(
        port,
        buffer_size
   );
     

}

#ifdef _WIN32
void PRECICE_MAINC2SOCKET_PLAIN_PORT_DESTROY_INSTANCE(long long *ptr){
#else
void precice_mainc2socket_plain_port_destroy_instance_(long long *ptr){
#endif
     precice::MainCxx2SocketPlainPort* c_ptr = (precice::MainCxx2SocketPlainPort*)*ptr;
     if(c_ptr!=NULL){
         delete c_ptr;
         c_ptr = NULL;
     }
}




#ifdef _WIN32
void PRECICE_MAINC2SOCKET_PLAIN_PORT_MAIN(long long* ref){
     
     ((precice::MainCxx2SocketPlainPort*)*ref)->main();
}
#else
void precice_mainc2socket_plain_port_main_(long long* ref){
     
     ((precice::MainCxx2SocketPlainPort*)*ref)->main();
}
#endif
}
