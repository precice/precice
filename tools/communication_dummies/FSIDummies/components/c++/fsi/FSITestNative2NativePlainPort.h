#ifndef FSI_FSITESTNATIVE2NATIVEPLAINPORT_H_
#define FSI_FSITESTNATIVE2NATIVEPLAINPORT_H_ 

#include "fsi/FSITest.h"
#include <jni.h> 
#include <iostream>

#ifdef __cplusplus
  extern "C" {
#endif


          
JNIEXPORT void JNICALL Java_fsi_FSITestNative2NativePlainPort_createInstance(JNIEnv *env, jobject obj);
JNIEXPORT void JNICALL Java_fsi_FSITestNative2NativePlainPort_destroyInstance(JNIEnv *env, jobject obj,jlong ref);
JNIEXPORT void JNICALL Java_fsi_FSITestNative2NativePlainPort_connect(JNIEnv *env, jobject obj,jlong ref,jlong destination);


#ifdef __cplusplus
  }
#endif

namespace fsi { 

     class FSITestNative2NativePlainPort;
}

class fsi::FSITestNative2NativePlainPort: public fsi::FSITest{
  private:
    fsi::FSITest* _destination;
  public:
    FSITestNative2NativePlainPort();
    ~FSITestNative2NativePlainPort();
    
    void connect(fsi::FSITest*);
    void test();  
    void testParallel();
   
};

#endif
