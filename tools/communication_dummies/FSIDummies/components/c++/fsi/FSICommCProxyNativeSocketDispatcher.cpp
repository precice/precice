
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
void FSI_FSICOMM_F2C_NSD_TRANSFERCOORDINATES(long long* ref,double* coord, int* coord_len){
#else
void fsi_fsicomm_f2c_nsd_transfercoordinates_(long long* ref,double* coord, int* coord_len){
#endif
    
     
     ((fsi::FSICommNativeSocketDispatcher*)*ref)->transferCoordinates(coord,*coord_len);
}
#ifdef _WIN32
void FSI_FSICOMM_F2C_NSD_TRANSFERDATA(long long* ref,double* data, int* data_len){
#else
void fsi_fsicomm_f2c_nsd_transferdata_(long long* ref,double* data, int* data_len){
#endif
    
     
     ((fsi::FSICommNativeSocketDispatcher*)*ref)->transferData(data,*data_len);
}


}
