#ifndef FSI_FSICOMMNATIVE2NATIVEPLAINPORT_H_
#define FSI_FSICOMMNATIVE2NATIVEPLAINPORT_H_ 

#include "fsi/FSIComm.h"
#include <jni.h> 
#include <iostream>

#ifdef __cplusplus
  extern "C" {
#endif


          
JNIEXPORT void JNICALL Java_fsi_FSICommNative2NativePlainPort_createInstance(JNIEnv *env, jobject obj);
JNIEXPORT void JNICALL Java_fsi_FSICommNative2NativePlainPort_destroyInstance(JNIEnv *env, jobject obj,jlong ref);
JNIEXPORT void JNICALL Java_fsi_FSICommNative2NativePlainPort_connect(JNIEnv *env, jobject obj,jlong ref,jlong destination);


#ifdef __cplusplus
  }
#endif

namespace fsi { 

     class FSICommNative2NativePlainPort;
}

class fsi::FSICommNative2NativePlainPort: public fsi::FSIComm{
  private:
    fsi::FSIComm* _destination;
  public:
    FSICommNative2NativePlainPort();
    ~FSICommNative2NativePlainPort();
    
    void connect(fsi::FSIComm*);
    void transferCoordinates(const int* coordId, const int coordId_len,const int* offsets, const int offsets_len,const std::string* hosts, const int hosts_len);  
    void transferCoordinatesParallel(const int* coordId, const int coordId_len,const int* offsets, const int offsets_len,const std::string* hosts, const int hosts_len);
   
    void transferData(const double* data, const int data_len);  
    void transferDataParallel(const double* data, const int data_len);
   
};

#endif
