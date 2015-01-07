#include "precice/ReceiverNative2NativePlainPort.h"
#include <assert.h>
#include "Component.h"

JNIEXPORT void JNICALL Java_precice_ReceiverNative2NativePlainPort_createInstance(JNIEnv *env, jobject obj){
  JavaVM* jvm;
  env->GetJavaVM(&jvm);
  
  precice::ReceiverNative2NativePlainPort *ref=new precice::ReceiverNative2NativePlainPort();
  
  jfieldID id =env->GetFieldID(env->GetObjectClass(obj), "_ref", "J");
  env->SetLongField(obj, id, (jlong)ref);
    
}

JNIEXPORT void JNICALL Java_precice_ReceiverNative2NativePlainPort_destroyInstance(JNIEnv *env, jobject obj,jlong ref){
  delete ((precice::ReceiverNative2NativePlainPort*)ref);
}

JNIEXPORT void JNICALL Java_precice_ReceiverNative2NativePlainPort_connect(JNIEnv *env, jobject obj,jlong ref,jlong destination){
  
  ((precice::ReceiverNative2NativePlainPort*)ref)->connect(dynamic_cast<precice::Receiver*>((Component*)destination));
}


precice::ReceiverNative2NativePlainPort::ReceiverNative2NativePlainPort():
     _destination(0){

}

precice::ReceiverNative2NativePlainPort::~ReceiverNative2NativePlainPort(){

}

void precice::ReceiverNative2NativePlainPort::connect(precice::Receiver* destination){
  _destination=destination;
}
void precice::ReceiverNative2NativePlainPort::receive(const double data,const int index,const int rank,int& tag){
     assert(_destination!=NULL);
     _destination->receive(data,index,rank,tag);
}
void precice::ReceiverNative2NativePlainPort::receiveParallel(const double data,const int index,const int rank,int& tag){
     assert(_destination!=NULL);
     _destination->receiveParallel(data,index,rank,tag);
}

