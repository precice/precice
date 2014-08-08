#include "fsi/FSIDataNative2NativePlainPort.h"
#include <assert.h>
#include "Component.h"

#ifdef JAVA
JNIEXPORT void JNICALL Java_fsi_FSIDataNative2NativePlainPort_createInstance(JNIEnv *env, jobject obj){
  JavaVM* jvm;
  env->GetJavaVM(&jvm);
  
  fsi::FSIDataNative2NativePlainPort *ref=new fsi::FSIDataNative2NativePlainPort();
  
  jfieldID id =env->GetFieldID(env->GetObjectClass(obj), "_ref", "J");
  env->SetLongField(obj, id, (jlong)ref);
    
}

JNIEXPORT void JNICALL Java_fsi_FSIDataNative2NativePlainPort_destroyInstance(JNIEnv *env, jobject obj,jlong ref){
  delete ((fsi::FSIDataNative2NativePlainPort*)ref);
}

JNIEXPORT void JNICALL Java_fsi_FSIDataNative2NativePlainPort_connect(JNIEnv *env, jobject obj,jlong ref,jlong destination){
  
  ((fsi::FSIDataNative2NativePlainPort*)ref)->connect(dynamic_cast<fsi::FSIData*>((Component*)destination));
}

#endif 

fsi::FSIDataNative2NativePlainPort::FSIDataNative2NativePlainPort():
     _destination(0){

}

fsi::FSIDataNative2NativePlainPort::~FSIDataNative2NativePlainPort(){

}

void fsi::FSIDataNative2NativePlainPort::connect(fsi::FSIData* destination){
  _destination=destination;
}
void fsi::FSIDataNative2NativePlainPort::transferData(const int* coordId, const int coordId_len,const double* data, const int data_len){
     assert(_destination!=NULL);
     _destination->transferData(coordId,coordId_len,data,data_len);
}
void fsi::FSIDataNative2NativePlainPort::transferDataParallel(const int* coordId, const int coordId_len,const double* data, const int data_len){
     assert(_destination!=NULL);
     _destination->transferDataParallel(coordId,coordId_len,data,data_len);
}
void fsi::FSIDataNative2NativePlainPort::dataAck(int& ack){
     assert(_destination!=NULL);
     _destination->dataAck(ack);
}
void fsi::FSIDataNative2NativePlainPort::dataAckParallel(int& ack){
     assert(_destination!=NULL);
     _destination->dataAckParallel(ack);
}

