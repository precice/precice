
#include <iostream>
#include <string.h>
#include "fsi/FSITestC2CxxSocketPlainPort.h"
#include "fsi/FSITestCxx2SocketPlainPort.h"
#include "fsi/FSITest.h"
extern "C" {

#ifdef _WIN32
void FSI_FSITESTC2SOCKET_PLAIN_PORT_CREATE_CLIENT_INSTANCE(long long* ptr,char* host,int& port,int& buffer_size){
#else
void fsi_fsitestc2socket_plain_port_create_client_instance_(long long* ptr,char* host,int& port,int& buffer_size){
#endif
   *ptr=(long long)new fsi::FSITestCxx2SocketPlainPort(
        host,
        port,
        buffer_size
   );
     

}

#ifdef _WIN32
void FSI_FSITESTC2SOCKET_PLAIN_PORT_CREATE_SERVER_INSTANCE(long long* ptr,int& port,int& buffer_size){
#else
void fsi_fsitestc2socket_plain_port_create_server_instance_(long long* ptr,int& port,int& buffer_size){
#endif
   *ptr=(long long)new fsi::FSITestCxx2SocketPlainPort(
        port,
        buffer_size
   );
     

}

#ifdef _WIN32
void FSI_FSITESTC2SOCKET_PLAIN_PORT_DESTROY_INSTANCE(long long *ptr){
#else
void fsi_fsitestc2socket_plain_port_destroy_instance_(long long *ptr){
#endif
     fsi::FSITestCxx2SocketPlainPort* c_ptr = (fsi::FSITestCxx2SocketPlainPort*)*ptr;
     if(c_ptr!=NULL){
         delete c_ptr;
         c_ptr = NULL;
     }
}




#ifdef _WIN32
void FSI_FSITESTC2SOCKET_PLAIN_PORT_TEST(long long* ref){
     
     ((fsi::FSITestCxx2SocketPlainPort*)*ref)->test();
}
#else
void fsi_fsitestc2socket_plain_port_test_(long long* ref){
     
     ((fsi::FSITestCxx2SocketPlainPort*)*ref)->test();
}
#endif
}
