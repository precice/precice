#include "fsi/FSIDataNativeDispatcher.h"
#include <algorithm>

#ifdef JAVA
JNIEXPORT void JNICALL Java_fsi_FSIDataNativeDispatcher_createInstance(JNIEnv *env, jobject obj){
  JavaVM* jvm;
  env->GetJavaVM(&jvm);
  
  fsi::FSIDataNativeDispatcher *ref=new fsi::FSIDataNativeDispatcher();
  
  jfieldID id =env->GetFieldID(env->GetObjectClass(obj), "_ref", "J");
  env->SetLongField(obj, id, (jlong)ref);
    
}

JNIEXPORT void JNICALL Java_fsi_FSIDataNativeDispatcher_destroyInstance(JNIEnv *env, jobject obj,jlong ref){
  delete ((fsi::FSIDataNativeDispatcher*)ref);
}

JNIEXPORT void JNICALL Java_fsi_FSIDataNativeDispatcher_connect(JNIEnv *env, jobject obj,jlong ref,jlong destination){
  ((fsi::FSIDataNativeDispatcher*)ref)->connect((fsi::FSIData*)destination);
}

JNIEXPORT void JNICALL Java_fsi_FSIDataNativeDispatcher_disconnect(JNIEnv *env, jobject obj,jlong ref,jlong destination){
  ((fsi::FSIDataNativeDispatcher*)ref)->disconnect((fsi::FSIData*)destination);
}
#endif

fsi::FSIDataNativeDispatcher::FSIDataNativeDispatcher(){

}

fsi::FSIDataNativeDispatcher::~FSIDataNativeDispatcher(){

}

void fsi::FSIDataNativeDispatcher::connect(fsi::FSIData* destination){
  if(std::find(_destinations.begin(), _destinations.end(), destination)==_destinations.end())
     _destinations.push_back(destination);
}

void fsi::FSIDataNativeDispatcher::disconnect(fsi::FSIData* destination){
  std::vector<fsi::FSIData*>::iterator iter=std::find(_destinations.begin(), _destinations.end(), destination);
  if(iter!=_destinations.end())
     _destinations.erase(iter);
}

bool fsi::FSIDataNativeDispatcher::isConnected() const{
  return !_destinations.empty();
}


void fsi::FSIDataNativeDispatcher::transferData(const int* coordId, const int coordId_len,const double* data, const int data_len){
    for(unsigned int i=0;i<_destinations.size();i++)
          _destinations[i]->transferData(coordId,coordId_len,data,data_len);
}

void fsi::FSIDataNativeDispatcher::transferDataParallel(const int* coordId, const int coordId_len,const double* data, const int data_len){
    for(unsigned int i=0;i<_destinations.size();i++)
          _destinations[i]->transferDataParallel(coordId,coordId_len,data,data_len);
}
void fsi::FSIDataNativeDispatcher::dataAck(int& ack){
    for(unsigned int i=0;i<_destinations.size();i++)
          _destinations[i]->dataAck(ack);
}

void fsi::FSIDataNativeDispatcher::dataAckParallel(int& ack){
    for(unsigned int i=0;i<_destinations.size();i++)
          _destinations[i]->dataAckParallel(ack);
}

