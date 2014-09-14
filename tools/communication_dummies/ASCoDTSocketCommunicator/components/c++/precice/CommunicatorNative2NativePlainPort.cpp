#include "precice/CommunicatorNative2NativePlainPort.h"
#include <assert.h>
#include "Component.h"

JNIEXPORT void JNICALL Java_precice_CommunicatorNative2NativePlainPort_createInstance(JNIEnv *env, jobject obj){
  JavaVM* jvm;
  env->GetJavaVM(&jvm);
  
  precice::CommunicatorNative2NativePlainPort *ref=new precice::CommunicatorNative2NativePlainPort();
  
  jfieldID id =env->GetFieldID(env->GetObjectClass(obj), "_ref", "J");
  env->SetLongField(obj, id, (jlong)ref);
    
}

JNIEXPORT void JNICALL Java_precice_CommunicatorNative2NativePlainPort_destroyInstance(JNIEnv *env, jobject obj,jlong ref){
  delete ((precice::CommunicatorNative2NativePlainPort*)ref);
}

JNIEXPORT void JNICALL Java_precice_CommunicatorNative2NativePlainPort_connect(JNIEnv *env, jobject obj,jlong ref,jlong destination){
  
  ((precice::CommunicatorNative2NativePlainPort*)ref)->connect(dynamic_cast<precice::Communicator*>((Component*)destination));
}


precice::CommunicatorNative2NativePlainPort::CommunicatorNative2NativePlainPort():
     _destination(0){

}

precice::CommunicatorNative2NativePlainPort::~CommunicatorNative2NativePlainPort(){

}

void precice::CommunicatorNative2NativePlainPort::connect(precice::Communicator* destination){
  _destination=destination;
}
void precice::CommunicatorNative2NativePlainPort::setData(const double* data, const int data_len){
     assert(_destination!=NULL);
     _destination->setData(data,data_len);
}
void precice::CommunicatorNative2NativePlainPort::setDataParallel(const double* data, const int data_len){
     assert(_destination!=NULL);
     _destination->setDataParallel(data,data_len);
}

