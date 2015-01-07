
#include <iostream>
#include <string.h>
#include "fsi/FSIDataC2CxxSocketPlainPort.h"
#include "fsi/FSIDataCxx2SocketPlainPort.h"
#include "fsi/FSIData.h"
extern "C" {

#ifdef _WIN32
void FSI_FSIDATAC2SOCKET_PLAIN_PORT_CREATE_CLIENT_INSTANCE(long long* ptr,char* host,int& port,int& buffer_size){
#else
void fsi_fsidatac2socket_plain_port_create_client_instance_(long long* ptr,char* host,int& port,int& buffer_size){
#endif
   *ptr=(long long)new fsi::FSIDataCxx2SocketPlainPort(
        host,
        port,
        buffer_size
   );
     

}

#ifdef _WIN32
void FSI_FSIDATAC2SOCKET_PLAIN_PORT_CREATE_SERVER_INSTANCE(long long* ptr,int& port,int& buffer_size){
#else
void fsi_fsidatac2socket_plain_port_create_server_instance_(long long* ptr,int& port,int& buffer_size){
#endif
   *ptr=(long long)new fsi::FSIDataCxx2SocketPlainPort(
        port,
        buffer_size
   );
     

}

#ifdef _WIN32
void FSI_FSIDATAC2SOCKET_PLAIN_PORT_DESTROY_INSTANCE(long long *ptr){
#else
void fsi_fsidatac2socket_plain_port_destroy_instance_(long long *ptr){
#endif
     fsi::FSIDataCxx2SocketPlainPort* c_ptr = (fsi::FSIDataCxx2SocketPlainPort*)*ptr;
     if(c_ptr!=NULL){
         delete c_ptr;
         c_ptr = NULL;
     }
}




#ifdef _WIN32
void FSI_FSIDATAC2SOCKET_PLAIN_PORT_TRANSFERDATA(long long* ref,int* coordId,int* coordId_len,double* data, int* data_len){
     
     ((fsi::FSIDataCxx2SocketPlainPort*)*ref)->transferData(coordId,*coordId_len,data,*data_len);
}
#else
void fsi_fsidatac2socket_plain_port_transferdata_(long long* ref,int* coordId,int* coordId_len,double* data, int* data_len){
     
     ((fsi::FSIDataCxx2SocketPlainPort*)*ref)->transferData(coordId,*coordId_len,data,*data_len);
}
#endif
#ifdef _WIN32
void FSI_FSIDATAC2SOCKET_PLAIN_PORT_DATAACK(long long* ref,int* ack){
     
     ((fsi::FSIDataCxx2SocketPlainPort*)*ref)->dataAck(*ack);
}
#else
void fsi_fsidatac2socket_plain_port_dataack_(long long* ref,int* ack){
     
     ((fsi::FSIDataCxx2SocketPlainPort*)*ref)->dataAck(*ack);
}
#endif
}
