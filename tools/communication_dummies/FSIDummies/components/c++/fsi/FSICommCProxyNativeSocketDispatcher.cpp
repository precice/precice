
#include "fsi/FSICommNativeSocketDispatcher.h"
#include "fsi/FSIComm.h"
#include <stdio.h>
#include <string.h>

extern "C" {

#ifdef _WIN32
void FSI_FSICOMM_F2C_NSD_CREATE_CLIENT_INSTANCE(long long* ptr,char* host,int& port,int& buffer_size){
#else
void fsi_fsicomm_f2c_nsd_create_client_instance_(long long* ptr,char* host,int& port,int& buffer_size){
#endif     
     *ptr=(long long)new fsi::FSICommNativeSocketDispatcher(
           host,
           port,
           buffer_size
     );
     

}




#ifdef _WIN32
void FSI_FSICOMM_F2C_NSD_DESTROY_INSTANCE(long long* ptr){

#else
void fsi_fsicomm_f2c_nsd_destroy_instance_(long long* ptr){
#endif
     fsi::FSICommNativeSocketDispatcher* dispatcher =
               ((fsi::FSICommNativeSocketDispatcher*)*ptr);
     if(dispatcher!=NULL){
          delete dispatcher;
          dispatcher=NULL;
     }
     

}

#ifdef _WIN32
void FSI_FSICOMM_F2C_NSD_TRANSFERCOORDINATES(long long* ref,int* coordId,int* coordId_len,int* offsets,int* offsets_len,char** hosts,int* hosts_len){
#else
void fsi_fsicomm_f2c_nsd_transfercoordinates_(long long* ref,int* coordId,int* coordId_len,int* offsets,int* offsets_len,char** hosts,int* hosts_len){
#endif
    
     std::string* hosts_str=new std::string[*hosts_len];
for(int i=0;i<*hosts_len;i++)
hosts_str[i]=hosts[i];

     ((fsi::FSICommNativeSocketDispatcher*)*ref)->transferCoordinates(coordId,*coordId_len,offsets,*offsets_len,hosts_str,*hosts_len);
}
#ifdef _WIN32
void FSI_FSICOMM_F2C_NSD_STARTDATATRANSFER(long long* ref){
#else
void fsi_fsicomm_f2c_nsd_startdatatransfer_(long long* ref){
#endif
    
     
     ((fsi::FSICommNativeSocketDispatcher*)*ref)->startDataTransfer();
}
#ifdef _WIN32
void FSI_FSICOMM_F2C_NSD_ENDDATATRANSFER(long long* ref,int* ack){
#else
void fsi_fsicomm_f2c_nsd_enddatatransfer_(long long* ref,int* ack){
#endif
    
     
     ((fsi::FSICommNativeSocketDispatcher*)*ref)->endDataTransfer(*ack);
}


}
