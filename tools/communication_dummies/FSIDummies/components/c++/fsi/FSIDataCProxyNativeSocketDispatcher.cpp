
#include "fsi/FSIDataNativeSocketDispatcher.h"
#include "fsi/FSIData.h"
#include <stdio.h>
#include <string.h>

extern "C" {

#ifdef _WIN32
void FSI_FSIDATA_F2C_NSD_CREATE_CLIENT_INSTANCE(long long* ptr,char* host,int& port,int& buffer_size){
#else
void fsi_fsidata_f2c_nsd_create_client_instance_(long long* ptr,char* host,int& port,int& buffer_size){
#endif     
     *ptr=(long long)new fsi::FSIDataNativeSocketDispatcher(
           host,
           port,
           buffer_size
     );
     

}




#ifdef _WIN32
void FSI_FSIDATA_F2C_NSD_DESTROY_INSTANCE(long long* ptr){

#else
void fsi_fsidata_f2c_nsd_destroy_instance_(long long* ptr){
#endif
     fsi::FSIDataNativeSocketDispatcher* dispatcher =
               ((fsi::FSIDataNativeSocketDispatcher*)*ptr);
     if(dispatcher!=NULL){
          delete dispatcher;
          dispatcher=NULL;
     }
     

}

#ifdef _WIN32
void FSI_FSIDATA_F2C_NSD_TRANSFERDATA(long long* ref,int* coordId,int* coordId_len,double* data, int* data_len){
#else
void fsi_fsidata_f2c_nsd_transferdata_(long long* ref,int* coordId,int* coordId_len,double* data, int* data_len){
#endif
    
     
     ((fsi::FSIDataNativeSocketDispatcher*)*ref)->transferData(coordId,*coordId_len,data,*data_len);
}
#ifdef _WIN32
void FSI_FSIDATA_F2C_NSD_DATAACK(long long* ref,int* ack){
#else
void fsi_fsidata_f2c_nsd_dataack_(long long* ref,int* ack){
#endif
    
     
     ((fsi::FSIDataNativeSocketDispatcher*)*ref)->dataAck(*ack);
}


}
