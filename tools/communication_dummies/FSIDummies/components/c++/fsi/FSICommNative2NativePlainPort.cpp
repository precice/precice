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
void fsi::FSICommNative2NativePlainPort::transferCoordinates(const int* coordId, const int coordId_len,const int* offsets, const int offsets_len,const std::string* hosts, const int hosts_len){
     assert(_destination!=NULL);
     _destination->transferCoordinates(coordId,coordId_len,offsets,offsets_len,hosts,hosts_len);
}
void fsi::FSICommNative2NativePlainPort::transferCoordinatesParallel(const int* coordId, const int coordId_len,const int* offsets, const int offsets_len,const std::string* hosts, const int hosts_len){
     assert(_destination!=NULL);
     _destination->transferCoordinatesParallel(coordId,coordId_len,offsets,offsets_len,hosts,hosts_len);
}
void fsi::FSICommNative2NativePlainPort::startDataTransfer(){
     assert(_destination!=NULL);
     _destination->startDataTransfer();
}
void fsi::FSICommNative2NativePlainPort::startDataTransferParallel(){
     assert(_destination!=NULL);
     _destination->startDataTransferParallel();
}
void fsi::FSICommNative2NativePlainPort::endDataTransfer(int& ack){
     assert(_destination!=NULL);
     _destination->endDataTransfer(ack);
}
void fsi::FSICommNative2NativePlainPort::endDataTransferParallel(int& ack){
     assert(_destination!=NULL);
     _destination->endDataTransferParallel(ack);
}

