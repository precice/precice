
#include "precice/ReceiverNativeSocketDispatcher.h"
#include "precice/Receiver.h"
#include <stdio.h>
#include <string.h>

extern "C" {

#ifdef _WIN32
void PRECICE_RECEIVER_F2C_NSD_CREATE_CLIENT_INSTANCE(long long* ptr,char* host,int& port,int& buffer_size){
#else
void precice_receiver_f2c_nsd_create_client_instance_(long long* ptr,char* host,int& port,int& buffer_size){
#endif     
     *ptr=(long long)new precice::ReceiverNativeSocketDispatcher(
           host,
           port,
           buffer_size
     );
     

}




#ifdef _WIN32
void PRECICE_RECEIVER_F2C_NSD_DESTROY_INSTANCE(long long* ptr){

#else
void precice_receiver_f2c_nsd_destroy_instance_(long long* ptr){
#endif
     precice::ReceiverNativeSocketDispatcher* dispatcher =
               ((precice::ReceiverNativeSocketDispatcher*)*ptr);
     if(dispatcher!=NULL){
          delete dispatcher;
          dispatcher=NULL;
     }
     

}

#ifdef _WIN32
void PRECICE_RECEIVER_F2C_NSD_RECEIVE(long long* ref,double* data,int* index,int* rank,int* tag){
#else
void precice_receiver_f2c_nsd_receive_(long long* ref,double* data,int* index,int* rank,int* tag){
#endif
    
     
     ((precice::ReceiverNativeSocketDispatcher*)*ref)->receive(*data,*index,*rank,*tag);
}


}
