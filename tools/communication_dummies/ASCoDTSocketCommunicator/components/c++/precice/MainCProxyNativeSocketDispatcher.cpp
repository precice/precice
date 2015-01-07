
#include "precice/MainNativeSocketDispatcher.h"
#include "precice/Main.h"
#include <stdio.h>
#include <string.h>

extern "C" {

#ifdef _WIN32
void PRECICE_MAIN_F2C_NSD_CREATE_CLIENT_INSTANCE(long long* ptr,char* host,int& port,int& buffer_size){
#else
void precice_main_f2c_nsd_create_client_instance_(long long* ptr,char* host,int& port,int& buffer_size){
#endif     
     *ptr=(long long)new precice::MainNativeSocketDispatcher(
           host,
           port,
           buffer_size
     );
     

}




#ifdef _WIN32
void PRECICE_MAIN_F2C_NSD_DESTROY_INSTANCE(long long* ptr){

#else
void precice_main_f2c_nsd_destroy_instance_(long long* ptr){
#endif
     precice::MainNativeSocketDispatcher* dispatcher =
               ((precice::MainNativeSocketDispatcher*)*ptr);
     if(dispatcher!=NULL){
          delete dispatcher;
          dispatcher=NULL;
     }
     

}

#ifdef _WIN32
void PRECICE_MAIN_F2C_NSD_MAIN(long long* ref){
#else
void precice_main_f2c_nsd_main_(long long* ref){
#endif
    
     
     ((precice::MainNativeSocketDispatcher*)*ref)->main();
}


}
