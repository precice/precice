
#include "fsi/FSITestNativeSocketDispatcher.h"
#include "fsi/FSITest.h"
#include <stdio.h>
#include <string.h>

extern "C" {

#ifdef _WIN32
void FSI_FSITEST_F2C_NSD_CREATE_CLIENT_INSTANCE(long long* ptr,char* host,int& port,int& buffer_size){
#else
void fsi_fsitest_f2c_nsd_create_client_instance_(long long* ptr,char* host,int& port,int& buffer_size){
#endif     
     *ptr=(long long)new fsi::FSITestNativeSocketDispatcher(
           host,
           port,
           buffer_size
     );
     

}




#ifdef _WIN32
void FSI_FSITEST_F2C_NSD_DESTROY_INSTANCE(long long* ptr){

#else
void fsi_fsitest_f2c_nsd_destroy_instance_(long long* ptr){
#endif
     fsi::FSITestNativeSocketDispatcher* dispatcher =
               ((fsi::FSITestNativeSocketDispatcher*)*ptr);
     if(dispatcher!=NULL){
          delete dispatcher;
          dispatcher=NULL;
     }
     

}

#ifdef _WIN32
void FSI_FSITEST_F2C_NSD_TEST(long long* ref){
#else
void fsi_fsitest_f2c_nsd_test_(long long* ref){
#endif
    
     
     ((fsi::FSITestNativeSocketDispatcher*)*ref)->test();
}


}
