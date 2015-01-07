#ifndef FSI_FSIDATANATIVE2NATIVEPLAINPORT_H_
#define FSI_FSIDATANATIVE2NATIVEPLAINPORT_H_ 

#include "fsi/FSIData.h"
#include <iostream>

#ifdef JAVA
#include <jni.h> 
#ifdef __cplusplus
  extern "C" {
#endif


          
JNIEXPORT void JNICALL Java_fsi_FSIDataNative2NativePlainPort_createInstance(JNIEnv *env, jobject obj);
JNIEXPORT void JNICALL Java_fsi_FSIDataNative2NativePlainPort_destroyInstance(JNIEnv *env, jobject obj,jlong ref);
JNIEXPORT void JNICALL Java_fsi_FSIDataNative2NativePlainPort_connect(JNIEnv *env, jobject obj,jlong ref,jlong destination);


#ifdef __cplusplus
  }
#endif
#endif

namespace fsi { 

     class FSIDataNative2NativePlainPort;
}

class fsi::FSIDataNative2NativePlainPort: public fsi::FSIData{
  private:
    fsi::FSIData* _destination;
  public:
    FSIDataNative2NativePlainPort();
    ~FSIDataNative2NativePlainPort();
    
    void connect(fsi::FSIData*);
    void transferData(const int* coordId, const int coordId_len,const double* data, const int data_len);  
    void transferDataParallel(const int* coordId, const int coordId_len,const double* data, const int data_len);
   
    void dataAck(int& ack);  
    void dataAckParallel(int& ack);
   
};

#endif
