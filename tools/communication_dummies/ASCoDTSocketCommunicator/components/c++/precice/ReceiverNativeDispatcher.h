#ifndef PRECICE_RECEIVERNATIVEDISPATCHER_H_
#define PRECICE_RECEIVERNATIVEDISPATCHER_H_ 

#include "precice/Receiver.h"
#include <iostream>
#include <vector>

namespace precice { 

     class ReceiverNativeDispatcher;
}

#include <jni.h> 

#ifdef __cplusplus
  extern "C" {
#endif


          
JNIEXPORT void JNICALL Java_precice_ReceiverNativeDispatcher_createInstance(JNIEnv *env, jobject obj);
JNIEXPORT void JNICALL Java_precice_ReceiverNativeDispatcher_destroyInstance(JNIEnv *env, jobject obj,jlong ref);
JNIEXPORT void JNICALL Java_precice_ReceiverNativeDispatcher_connect(JNIEnv *env, jobject obj,jlong ref,jlong port);
JNIEXPORT void JNICALL Java_precice_ReceiverNativeDispatcher_disconnect(JNIEnv *env, jobject obj,jlong ref,jlong port);


#ifdef __cplusplus
  }
#endif

class precice::ReceiverNativeDispatcher: public precice::Receiver{
  protected:
    std::vector<precice::Receiver*> _destinations;
  public:
    ReceiverNativeDispatcher();
    virtual ~ReceiverNativeDispatcher();
    
    void connect(precice::Receiver* ref);
    void disconnect(precice::Receiver* ref);
    bool isConnected() const;
    void receive(const double data,const int index,const int rank,int& tag);  
    void receiveParallel(const double data,const int index,const int rank,int& tag);
   
};

#endif
