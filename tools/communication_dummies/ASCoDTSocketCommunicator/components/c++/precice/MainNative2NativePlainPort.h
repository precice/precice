#ifndef PRECICE_MAINNATIVE2NATIVEPLAINPORT_H_
#define PRECICE_MAINNATIVE2NATIVEPLAINPORT_H_ 

#include "precice/Main.h"
#include <jni.h> 
#include <iostream>

#ifdef __cplusplus
  extern "C" {
#endif


          
JNIEXPORT void JNICALL Java_precice_MainNative2NativePlainPort_createInstance(JNIEnv *env, jobject obj);
JNIEXPORT void JNICALL Java_precice_MainNative2NativePlainPort_destroyInstance(JNIEnv *env, jobject obj,jlong ref);
JNIEXPORT void JNICALL Java_precice_MainNative2NativePlainPort_connect(JNIEnv *env, jobject obj,jlong ref,jlong destination);


#ifdef __cplusplus
  }
#endif

namespace precice { 

     class MainNative2NativePlainPort;
}

class precice::MainNative2NativePlainPort: public precice::Main{
  private:
    precice::Main* _destination;
  public:
    MainNative2NativePlainPort();
    ~MainNative2NativePlainPort();
    
    void connect(precice::Main*);
    void main();  
    void mainParallel();
   
};

#endif
