#ifndef PRECICE_INITIALIZERNATIVEDISPATCHER_H_
#define PRECICE_INITIALIZERNATIVEDISPATCHER_H_ 

#include "precice/Initializer.h"
#include <iostream>
#include <vector>

namespace precice { 

     class InitializerNativeDispatcher;
}

#include <jni.h> 

#ifdef __cplusplus
  extern "C" {
#endif


          
JNIEXPORT void JNICALL Java_precice_InitializerNativeDispatcher_createInstance(JNIEnv *env, jobject obj);
JNIEXPORT void JNICALL Java_precice_InitializerNativeDispatcher_destroyInstance(JNIEnv *env, jobject obj,jlong ref);
JNIEXPORT void JNICALL Java_precice_InitializerNativeDispatcher_connect(JNIEnv *env, jobject obj,jlong ref,jlong port);
JNIEXPORT void JNICALL Java_precice_InitializerNativeDispatcher_disconnect(JNIEnv *env, jobject obj,jlong ref,jlong port);


#ifdef __cplusplus
  }
#endif

class precice::InitializerNativeDispatcher: public precice::Initializer{
  protected:
    std::vector<precice::Initializer*> _destinations;
  public:
    InitializerNativeDispatcher();
    virtual ~InitializerNativeDispatcher();
    
    void connect(precice::Initializer* ref);
    void disconnect(precice::Initializer* ref);
    bool isConnected() const;
    void acknowledge(const int identifier,int& tag);  
    void acknowledgeParallel(const int identifier,int& tag);
   
    void initialize(const std::string* addresses, const int addresses_len,const int* vertexes, const int vertexes_len);  
    void initializeParallel(const std::string* addresses, const int addresses_len,const int* vertexes, const int vertexes_len);
   
};

#endif
