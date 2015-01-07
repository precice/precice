#include "precice/MainNative2NativePlainPort.h"
#include <assert.h>
#include "Component.h"

JNIEXPORT void JNICALL Java_precice_MainNative2NativePlainPort_createInstance(JNIEnv *env, jobject obj){
  JavaVM* jvm;
  env->GetJavaVM(&jvm);
  
  precice::MainNative2NativePlainPort *ref=new precice::MainNative2NativePlainPort();
  
  jfieldID id =env->GetFieldID(env->GetObjectClass(obj), "_ref", "J");
  env->SetLongField(obj, id, (jlong)ref);
    
}

JNIEXPORT void JNICALL Java_precice_MainNative2NativePlainPort_destroyInstance(JNIEnv *env, jobject obj,jlong ref){
  delete ((precice::MainNative2NativePlainPort*)ref);
}

JNIEXPORT void JNICALL Java_precice_MainNative2NativePlainPort_connect(JNIEnv *env, jobject obj,jlong ref,jlong destination){
  
  ((precice::MainNative2NativePlainPort*)ref)->connect(dynamic_cast<precice::Main*>((Component*)destination));
}


precice::MainNative2NativePlainPort::MainNative2NativePlainPort():
     _destination(0){

}

precice::MainNative2NativePlainPort::~MainNative2NativePlainPort(){

}

void precice::MainNative2NativePlainPort::connect(precice::Main* destination){
  _destination=destination;
}
void precice::MainNative2NativePlainPort::main(){
     assert(_destination!=NULL);
     _destination->main();
}
void precice::MainNative2NativePlainPort::mainParallel(){
     assert(_destination!=NULL);
     _destination->mainParallel();
}

