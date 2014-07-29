
#include <iostream>
#include <string.h>
#include "fsi/FSICommC2CxxSocketPlainPort.h"
#include "fsi/FSICommCxx2SocketPlainPort.h"
#include "fsi/FSIComm.h"
extern "C" {

#ifdef _WIN32
void FSI_FSICOMMC2SOCKET_PLAIN_PORT_CREATE_CLIENT_INSTANCE(long long* ptr,char* host,int& port,int& buffer_size){
#else
void fsi_fsicommc2socket_plain_port_create_client_instance_(long long* ptr,char* host,int& port,int& buffer_size){
#endif
   *ptr=(long long)new fsi::FSICommCxx2SocketPlainPort(
        host,
        port,
        buffer_size
   );
     

}

#ifdef _WIN32
void FSI_FSICOMMC2SOCKET_PLAIN_PORT_CREATE_SERVER_INSTANCE(long long* ptr,int& port,int& buffer_size){
#else
void fsi_fsicommc2socket_plain_port_create_server_instance_(long long* ptr,int& port,int& buffer_size){
#endif
   *ptr=(long long)new fsi::FSICommCxx2SocketPlainPort(
        port,
        buffer_size
   );
     

}

#ifdef _WIN32
void FSI_FSICOMMC2SOCKET_PLAIN_PORT_DESTROY_INSTANCE(long long *ptr){
#else
void fsi_fsicommc2socket_plain_port_destroy_instance_(long long *ptr){
#endif
     fsi::FSICommCxx2SocketPlainPort* c_ptr = (fsi::FSICommCxx2SocketPlainPort*)*ptr;
     if(c_ptr!=NULL){
         delete c_ptr;
         c_ptr = NULL;
     }
}




#ifdef _WIN32
void FSI_FSICOMMC2SOCKET_PLAIN_PORT_TRANSFERCOORDINATES(long long* ref,double* coord, int* coord_len){
     
     ((fsi::FSICommCxx2SocketPlainPort*)*ref)->transferCoordinates(coord,*coord_len);
}
#else
void fsi_fsicommc2socket_plain_port_transfercoordinates_(long long* ref,double* coord, int* coord_len){
     
     ((fsi::FSICommCxx2SocketPlainPort*)*ref)->transferCoordinates(coord,*coord_len);
}
#endif
#ifdef _WIN32
void FSI_FSICOMMC2SOCKET_PLAIN_PORT_TRANSFERDATA(long long* ref,double* data, int* data_len){
     
     ((fsi::FSICommCxx2SocketPlainPort*)*ref)->transferData(data,*data_len);
}
#else
void fsi_fsicommc2socket_plain_port_transferdata_(long long* ref,double* data, int* data_len){
     
     ((fsi::FSICommCxx2SocketPlainPort*)*ref)->transferData(data,*data_len);
}
#endif
}
