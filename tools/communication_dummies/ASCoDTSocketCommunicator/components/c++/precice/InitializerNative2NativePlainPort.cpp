#include "precice/InitializerNative2NativePlainPort.h"
#include <assert.h>
#include "Component.h"

JNIEXPORT void JNICALL Java_precice_InitializerNative2NativePlainPort_createInstance(JNIEnv *env, jobject obj){
  JavaVM* jvm;
  env->GetJavaVM(&jvm);
  
  precice::InitializerNative2NativePlainPort *ref=new precice::InitializerNative2NativePlainPort();
  
  jfieldID id =env->GetFieldID(env->GetObjectClass(obj), "_ref", "J");
  env->SetLongField(obj, id, (jlong)ref);
    
}

JNIEXPORT void JNICALL Java_precice_InitializerNative2NativePlainPort_destroyInstance(JNIEnv *env, jobject obj,jlong ref){
  delete ((precice::InitializerNative2NativePlainPort*)ref);
}

JNIEXPORT void JNICALL Java_precice_InitializerNative2NativePlainPort_connect(JNIEnv *env, jobject obj,jlong ref,jlong destination){
  
  ((precice::InitializerNative2NativePlainPort*)ref)->connect(dynamic_cast<precice::Initializer*>((Component*)destination));
}


precice::InitializerNative2NativePlainPort::InitializerNative2NativePlainPort():
     _destination(0){

}

precice::InitializerNative2NativePlainPort::~InitializerNative2NativePlainPort(){

}

void precice::InitializerNative2NativePlainPort::connect(precice::Initializer* destination){
  _destination=destination;
}
void precice::InitializerNative2NativePlainPort::acknowledge(const int identifier,int& tag){
     assert(_destination!=NULL);
     _destination->acknowledge(identifier,tag);
}
void precice::InitializerNative2NativePlainPort::acknowledgeParallel(const int identifier,int& tag){
     assert(_destination!=NULL);
     _destination->acknowledgeParallel(identifier,tag);
}
void precice::InitializerNative2NativePlainPort::initialize(const std::string* addresses, const int addresses_len,const int* vertexes, const int vertexes_len){
     assert(_destination!=NULL);
     _destination->initialize(addresses,addresses_len,vertexes,vertexes_len);
}
void precice::InitializerNative2NativePlainPort::initializeParallel(const std::string* addresses, const int addresses_len,const int* vertexes, const int vertexes_len){
     assert(_destination!=NULL);
     _destination->initializeParallel(addresses,addresses_len,vertexes,vertexes_len);
}

