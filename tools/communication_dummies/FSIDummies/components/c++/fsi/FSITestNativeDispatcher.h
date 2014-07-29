#ifndef FSI_FSITESTNATIVEDISPATCHER_H_
#define FSI_FSITESTNATIVEDISPATCHER_H_ 

#include "fsi/FSITest.h"
#include <iostream>
#include <vector>

namespace fsi { 

     class FSITestNativeDispatcher;
}

#include <jni.h> 

#ifdef __cplusplus
  extern "C" {
#endif


          
JNIEXPORT void JNICALL Java_fsi_FSITestNativeDispatcher_createInstance(JNIEnv *env, jobject obj);
JNIEXPORT void JNICALL Java_fsi_FSITestNativeDispatcher_destroyInstance(JNIEnv *env, jobject obj,jlong ref);
JNIEXPORT void JNICALL Java_fsi_FSITestNativeDispatcher_connect(JNIEnv *env, jobject obj,jlong ref,jlong port);
JNIEXPORT void JNICALL Java_fsi_FSITestNativeDispatcher_disconnect(JNIEnv *env, jobject obj,jlong ref,jlong port);


#ifdef __cplusplus
  }
#endif

class fsi::FSITestNativeDispatcher: public fsi::FSITest{
  protected:
    std::vector<fsi::FSITest*> _destinations;
  public:
    FSITestNativeDispatcher();
    virtual ~FSITestNativeDispatcher();
    
    void connect(fsi::FSITest* ref);
    void disconnect(fsi::FSITest* ref);
    bool isConnected() const;
    void test();  
    void testParallel();
   
};

#endif
