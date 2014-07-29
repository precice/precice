#ifndef FSI_FSICOMMNATIVEDISPATCHER_H_
#define FSI_FSICOMMNATIVEDISPATCHER_H_ 

#include "fsi/FSIComm.h"
#include <iostream>
#include <vector>

namespace fsi { 

     class FSICommNativeDispatcher;
}

#include <jni.h> 

#ifdef __cplusplus
  extern "C" {
#endif


          
JNIEXPORT void JNICALL Java_fsi_FSICommNativeDispatcher_createInstance(JNIEnv *env, jobject obj);
JNIEXPORT void JNICALL Java_fsi_FSICommNativeDispatcher_destroyInstance(JNIEnv *env, jobject obj,jlong ref);
JNIEXPORT void JNICALL Java_fsi_FSICommNativeDispatcher_connect(JNIEnv *env, jobject obj,jlong ref,jlong port);
JNIEXPORT void JNICALL Java_fsi_FSICommNativeDispatcher_disconnect(JNIEnv *env, jobject obj,jlong ref,jlong port);


#ifdef __cplusplus
  }
#endif

class fsi::FSICommNativeDispatcher: public fsi::FSIComm{
  protected:
    std::vector<fsi::FSIComm*> _destinations;
  public:
    FSICommNativeDispatcher();
    virtual ~FSICommNativeDispatcher();
    
    void connect(fsi::FSIComm* ref);
    void disconnect(fsi::FSIComm* ref);
    bool isConnected() const;
    void transferCoordinates(const double* coord, const int coord_len);  
    void transferCoordinatesParallel(const double* coord, const int coord_len);
   
    void transferData(const double* data, const int data_len);  
    void transferDataParallel(const double* data, const int data_len);
   
};

#endif
