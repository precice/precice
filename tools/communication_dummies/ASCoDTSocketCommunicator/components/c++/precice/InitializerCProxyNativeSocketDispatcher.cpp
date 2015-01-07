
#include "precice/InitializerNativeSocketDispatcher.h"
#include "precice/Initializer.h"
#include <stdio.h>
#include <string.h>

extern "C" {

#ifdef _WIN32
void PRECICE_INITIALIZER_F2C_NSD_CREATE_CLIENT_INSTANCE(long long* ptr,char* host,int& port,int& buffer_size){
#else
void precice_initializer_f2c_nsd_create_client_instance_(long long* ptr,char* host,int& port,int& buffer_size){
#endif     
     *ptr=(long long)new precice::InitializerNativeSocketDispatcher(
           host,
           port,
           buffer_size
     );
     

}




#ifdef _WIN32
void PRECICE_INITIALIZER_F2C_NSD_DESTROY_INSTANCE(long long* ptr){

#else
void precice_initializer_f2c_nsd_destroy_instance_(long long* ptr){
#endif
     precice::InitializerNativeSocketDispatcher* dispatcher =
               ((precice::InitializerNativeSocketDispatcher*)*ptr);
     if(dispatcher!=NULL){
          delete dispatcher;
          dispatcher=NULL;
     }
     

}

#ifdef _WIN32
void PRECICE_INITIALIZER_F2C_NSD_ACKNOWLEDGE(long long* ref,int* identifier,int* tag){
#else
void precice_initializer_f2c_nsd_acknowledge_(long long* ref,int* identifier,int* tag){
#endif
    
     
     ((precice::InitializerNativeSocketDispatcher*)*ref)->acknowledge(*identifier,*tag);
}
#ifdef _WIN32
void PRECICE_INITIALIZER_F2C_NSD_INITIALIZE(long long* ref,char** addresses,int* addresses_len,int* vertexes,int* vertexes_len){
#else
void precice_initializer_f2c_nsd_initialize_(long long* ref,char** addresses,int* addresses_len,int* vertexes,int* vertexes_len){
#endif
    
     std::string* addresses_str=new std::string[*addresses_len];
for(int i=0;i<*addresses_len;i++)
addresses_str[i]=addresses[i];

     ((precice::InitializerNativeSocketDispatcher*)*ref)->initialize(addresses_str,*addresses_len,vertexes,*vertexes_len);
}


}
