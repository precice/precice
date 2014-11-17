#ifndef PRECICE_RECEIVERNATIVE2NATIVEPLAINPORT_H_
#define PRECICE_RECEIVERNATIVE2NATIVEPLAINPORT_H_ 

#include "precice/Receiver.h"
#include <jni.h> 
#include <iostream>

#ifdef __cplusplus
  extern "C" {
#endif


          
JNIEXPORT void JNICALL Java_precice_ReceiverNative2NativePlainPort_createInstance(JNIEnv *env, jobject obj);
JNIEXPORT void JNICALL Java_precice_ReceiverNative2NativePlainPort_destroyInstance(JNIEnv *env, jobject obj,jlong ref);
JNIEXPORT void JNICALL Java_precice_ReceiverNative2NativePlainPort_connect(JNIEnv *env, jobject obj,jlong ref,jlong destination);


#ifdef __cplusplus
  }
#endif

namespace precice { 

     class ReceiverNative2NativePlainPort;
}

class precice::ReceiverNative2NativePlainPort: public precice::Receiver{
  private:
    precice::Receiver* _destination;
  public:
    ReceiverNative2NativePlainPort();
    ~ReceiverNative2NativePlainPort();
    
    void connect(precice::Receiver*);
    void receive(const double data,const int index,const int rank,int& tag);  
    void receiveParallel(const double data,const int index,const int rank,int& tag);
   
};

#endif
