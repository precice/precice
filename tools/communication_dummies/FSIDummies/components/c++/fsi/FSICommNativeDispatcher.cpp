#include "fsi/FSICommNativeDispatcher.h"
#include <algorithm>

JNIEXPORT void JNICALL Java_fsi_FSICommNativeDispatcher_createInstance(JNIEnv *env, jobject obj){
  JavaVM* jvm;
  env->GetJavaVM(&jvm);
  
  fsi::FSICommNativeDispatcher *ref=new fsi::FSICommNativeDispatcher();
  
  jfieldID id =env->GetFieldID(env->GetObjectClass(obj), "_ref", "J");
  env->SetLongField(obj, id, (jlong)ref);
    
}

JNIEXPORT void JNICALL Java_fsi_FSICommNativeDispatcher_destroyInstance(JNIEnv *env, jobject obj,jlong ref){
  delete ((fsi::FSICommNativeDispatcher*)ref);
}

JNIEXPORT void JNICALL Java_fsi_FSICommNativeDispatcher_connect(JNIEnv *env, jobject obj,jlong ref,jlong destination){
  ((fsi::FSICommNativeDispatcher*)ref)->connect((fsi::FSIComm*)destination);
}

JNIEXPORT void JNICALL Java_fsi_FSICommNativeDispatcher_disconnect(JNIEnv *env, jobject obj,jlong ref,jlong destination){
  ((fsi::FSICommNativeDispatcher*)ref)->disconnect((fsi::FSIComm*)destination);
}


fsi::FSICommNativeDispatcher::FSICommNativeDispatcher(){

}

fsi::FSICommNativeDispatcher::~FSICommNativeDispatcher(){

}

void fsi::FSICommNativeDispatcher::connect(fsi::FSIComm* destination){
  if(std::find(_destinations.begin(), _destinations.end(), destination)==_destinations.end())
     _destinations.push_back(destination);
}

void fsi::FSICommNativeDispatcher::disconnect(fsi::FSIComm* destination){
  std::vector<fsi::FSIComm*>::iterator iter=std::find(_destinations.begin(), _destinations.end(), destination);
  if(iter!=_destinations.end())
     _destinations.erase(iter);
}

bool fsi::FSICommNativeDispatcher::isConnected() const{
  return !_destinations.empty();
}


void fsi::FSICommNativeDispatcher::transferCoordinates(const double* coord, const int coord_len){
    for(unsigned int i=0;i<_destinations.size();i++)
          _destinations[i]->transferCoordinates(coord,coord_len);
}

void fsi::FSICommNativeDispatcher::transferCoordinatesParallel(const double* coord, const int coord_len){
    for(unsigned int i=0;i<_destinations.size();i++)
          _destinations[i]->transferCoordinatesParallel(coord,coord_len);
}
void fsi::FSICommNativeDispatcher::transferData(const double* data, const int data_len){
    for(unsigned int i=0;i<_destinations.size();i++)
          _destinations[i]->transferData(data,data_len);
}

void fsi::FSICommNativeDispatcher::transferDataParallel(const double* data, const int data_len){
    for(unsigned int i=0;i<_destinations.size();i++)
          _destinations[i]->transferDataParallel(data,data_len);
}

