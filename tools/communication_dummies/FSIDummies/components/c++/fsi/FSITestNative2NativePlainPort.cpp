#include "fsi/FSITestNative2NativePlainPort.h"
#include <assert.h>
#include "Component.h"

JNIEXPORT void JNICALL Java_fsi_FSITestNative2NativePlainPort_createInstance(JNIEnv *env, jobject obj){
  JavaVM* jvm;
  env->GetJavaVM(&jvm);
  
  fsi::FSITestNative2NativePlainPort *ref=new fsi::FSITestNative2NativePlainPort();
  
  jfieldID id =env->GetFieldID(env->GetObjectClass(obj), "_ref", "J");
  env->SetLongField(obj, id, (jlong)ref);
    
}

JNIEXPORT void JNICALL Java_fsi_FSITestNative2NativePlainPort_destroyInstance(JNIEnv *env, jobject obj,jlong ref){
  delete ((fsi::FSITestNative2NativePlainPort*)ref);
}

JNIEXPORT void JNICALL Java_fsi_FSITestNative2NativePlainPort_connect(JNIEnv *env, jobject obj,jlong ref,jlong destination){
  
  ((fsi::FSITestNative2NativePlainPort*)ref)->connect(dynamic_cast<fsi::FSITest*>((Component*)destination));
}


fsi::FSITestNative2NativePlainPort::FSITestNative2NativePlainPort():
     _destination(0){

}

fsi::FSITestNative2NativePlainPort::~FSITestNative2NativePlainPort(){

}

void fsi::FSITestNative2NativePlainPort::connect(fsi::FSITest* destination){
  _destination=destination;
}
void fsi::FSITestNative2NativePlainPort::test(){
     assert(_destination!=NULL);
     _destination->test();
}
void fsi::FSITestNative2NativePlainPort::testParallel(){
     assert(_destination!=NULL);
     _destination->testParallel();
}

