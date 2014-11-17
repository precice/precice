#ifndef PRECICE_MAINNATIVEDISPATCHER_H_
#define PRECICE_MAINNATIVEDISPATCHER_H_ 

#include "precice/Main.h"
#include <iostream>
#include <vector>

namespace precice { 

     class MainNativeDispatcher;
}

#include <jni.h> 

#ifdef __cplusplus
  extern "C" {
#endif


          
JNIEXPORT void JNICALL Java_precice_MainNativeDispatcher_createInstance(JNIEnv *env, jobject obj);
JNIEXPORT void JNICALL Java_precice_MainNativeDispatcher_destroyInstance(JNIEnv *env, jobject obj,jlong ref);
JNIEXPORT void JNICALL Java_precice_MainNativeDispatcher_connect(JNIEnv *env, jobject obj,jlong ref,jlong port);
JNIEXPORT void JNICALL Java_precice_MainNativeDispatcher_disconnect(JNIEnv *env, jobject obj,jlong ref,jlong port);


#ifdef __cplusplus
  }
#endif

class precice::MainNativeDispatcher: public precice::Main{
  protected:
    std::vector<precice::Main*> _destinations;
  public:
    MainNativeDispatcher();
    virtual ~MainNativeDispatcher();
    
    void connect(precice::Main* ref);
    void disconnect(precice::Main* ref);
    bool isConnected() const;
    void main();  
    void mainParallel();
   
};

#endif
