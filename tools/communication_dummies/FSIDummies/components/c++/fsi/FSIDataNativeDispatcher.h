#ifndef FSI_FSIDATANATIVEDISPATCHER_H_
#define FSI_FSIDATANATIVEDISPATCHER_H_ 

#include "fsi/FSIData.h"
#include <iostream>
#include <vector>

namespace fsi { 

     class FSIDataNativeDispatcher;
}

#ifdef JAVA
#include <jni.h> 

#ifdef __cplusplus
  extern "C" {
#endif


          
JNIEXPORT void JNICALL Java_fsi_FSIDataNativeDispatcher_createInstance(JNIEnv *env, jobject obj);
JNIEXPORT void JNICALL Java_fsi_FSIDataNativeDispatcher_destroyInstance(JNIEnv *env, jobject obj,jlong ref);
JNIEXPORT void JNICALL Java_fsi_FSIDataNativeDispatcher_connect(JNIEnv *env, jobject obj,jlong ref,jlong port);
JNIEXPORT void JNICALL Java_fsi_FSIDataNativeDispatcher_disconnect(JNIEnv *env, jobject obj,jlong ref,jlong port);


#ifdef __cplusplus
  }
#endif

#endif

class fsi::FSIDataNativeDispatcher: public fsi::FSIData{
  protected:
    std::vector<fsi::FSIData*> _destinations;
  public:
    FSIDataNativeDispatcher();
    virtual ~FSIDataNativeDispatcher();
    
    void connect(fsi::FSIData* ref);
    void disconnect(fsi::FSIData* ref);
    bool isConnected() const;
    void transferData(const int* coordId, const int coordId_len,const double* data, const int data_len);  
    void transferDataParallel(const int* coordId, const int coordId_len,const double* data, const int data_len);
   
    void dataAck(int& ack);  
    void dataAckParallel(int& ack);
   
};

#endif
