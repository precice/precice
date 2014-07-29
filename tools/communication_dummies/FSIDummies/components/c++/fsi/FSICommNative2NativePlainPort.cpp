#include "fsi/FSICommNative2NativePlainPort.h"
#include <assert.h>
#include "Component.h"

JNIEXPORT void JNICALL Java_fsi_FSICommNative2NativePlainPort_createInstance(JNIEnv *env, jobject obj){
  JavaVM* jvm;
  env->GetJavaVM(&jvm);
  
  fsi::FSICommNative2NativePlainPort *ref=new fsi::FSICommNative2NativePlainPort();
  
  jfieldID id =env->GetFieldID(env->GetObjectClass(obj), "_ref", "J");
  env->SetLongField(obj, id, (jlong)ref);
    
}

JNIEXPORT void JNICALL Java_fsi_FSICommNative2NativePlainPort_destroyInstance(JNIEnv *env, jobject obj,jlong ref){
  delete ((fsi::FSICommNative2NativePlainPort*)ref);
}

JNIEXPORT void JNICALL Java_fsi_FSICommNative2NativePlainPort_connect(JNIEnv *env, jobject obj,jlong ref,jlong destination){
  
  ((fsi::FSICommNative2NativePlainPort*)ref)->connect(dynamic_cast<fsi::FSIComm*>((Component*)destination));
}


fsi::FSICommNative2NativePlainPort::FSICommNative2NativePlainPort():
     _destination(0){

}

fsi::FSICommNative2NativePlainPort::~FSICommNative2NativePlainPort(){

}

void fsi::FSICommNative2NativePlainPort::connect(fsi::FSIComm* destination){
  _destination=destination;
}
void fsi::FSICommNative2NativePlainPort::transferCoordinates(const double* coord, const int coord_len){
     assert(_destination!=NULL);
     _destination->transferCoordinates(coord,coord_len);
}
void fsi::FSICommNative2NativePlainPort::transferCoordinatesParallel(const double* coord, const int coord_len){
     assert(_destination!=NULL);
     _destination->transferCoordinatesParallel(coord,coord_len);
}
void fsi::FSICommNative2NativePlainPort::transferData(const double* data, const int data_len){
     assert(_destination!=NULL);
     _destination->transferData(data,data_len);
}
void fsi::FSICommNative2NativePlainPort::transferDataParallel(const double* data, const int data_len){
     assert(_destination!=NULL);
     _destination->transferDataParallel(data,data_len);
}

