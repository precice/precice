#ifndef PRECICE_COMMUNICATORNATIVE2NATIVEPLAINPORT_H_
#define PRECICE_COMMUNICATORNATIVE2NATIVEPLAINPORT_H_ 

#include "precice/Communicator.h"
#include <jni.h> 
#include <iostream>

#ifdef __cplusplus
  extern "C" {
#endif


          
JNIEXPORT void JNICALL Java_precice_CommunicatorNative2NativePlainPort_createInstance(JNIEnv *env, jobject obj);
JNIEXPORT void JNICALL Java_precice_CommunicatorNative2NativePlainPort_destroyInstance(JNIEnv *env, jobject obj,jlong ref);
JNIEXPORT void JNICALL Java_precice_CommunicatorNative2NativePlainPort_connect(JNIEnv *env, jobject obj,jlong ref,jlong destination);


#ifdef __cplusplus
  }
#endif

namespace precice { 

     class CommunicatorNative2NativePlainPort;
}

class precice::CommunicatorNative2NativePlainPort: public precice::Communicator{
  private:
    precice::Communicator* _destination;
  public:
    CommunicatorNative2NativePlainPort();
    ~CommunicatorNative2NativePlainPort();
    
    void connect(precice::Communicator*);
    void setData(const double data,const int index,const int rank,int& tag);  
    void setDataParallel(const double data,const int index,const int rank,int& tag);
   
};

#endif
