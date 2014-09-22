#ifndef PRECICE_COMMUNICATORNATIVEDISPATCHER_H_
#define PRECICE_COMMUNICATORNATIVEDISPATCHER_H_ 

#include "precice/Communicator.h"
#include <iostream>
#include <vector>

namespace precice { 

     class CommunicatorNativeDispatcher;
}

#include <jni.h> 

#ifdef __cplusplus
  extern "C" {
#endif


          
JNIEXPORT void JNICALL Java_precice_CommunicatorNativeDispatcher_createInstance(JNIEnv *env, jobject obj);
JNIEXPORT void JNICALL Java_precice_CommunicatorNativeDispatcher_destroyInstance(JNIEnv *env, jobject obj,jlong ref);
JNIEXPORT void JNICALL Java_precice_CommunicatorNativeDispatcher_connect(JNIEnv *env, jobject obj,jlong ref,jlong port);
JNIEXPORT void JNICALL Java_precice_CommunicatorNativeDispatcher_disconnect(JNIEnv *env, jobject obj,jlong ref,jlong port);


#ifdef __cplusplus
  }
#endif

class precice::CommunicatorNativeDispatcher: public precice::Communicator{
  protected:
    std::vector<precice::Communicator*> _destinations;
  public:
    CommunicatorNativeDispatcher();
    virtual ~CommunicatorNativeDispatcher();
    
    void connect(precice::Communicator* ref);
    void disconnect(precice::Communicator* ref);
    bool isConnected() const;
    void setData(const double data,const int index,const int rank,int& tag);  
    void setDataParallel(const double data,const int index,const int rank,int& tag);
   
};

#endif
