#include "precice/ReceiverNativeDispatcher.h"
#include <algorithm>

JNIEXPORT void JNICALL Java_precice_ReceiverNativeDispatcher_createInstance(JNIEnv *env, jobject obj){
  JavaVM* jvm;
  env->GetJavaVM(&jvm);
  
  precice::ReceiverNativeDispatcher *ref=new precice::ReceiverNativeDispatcher();
  
  jfieldID id =env->GetFieldID(env->GetObjectClass(obj), "_ref", "J");
  env->SetLongField(obj, id, (jlong)ref);
    
}

JNIEXPORT void JNICALL Java_precice_ReceiverNativeDispatcher_destroyInstance(JNIEnv *env, jobject obj,jlong ref){
  delete ((precice::ReceiverNativeDispatcher*)ref);
}

JNIEXPORT void JNICALL Java_precice_ReceiverNativeDispatcher_connect(JNIEnv *env, jobject obj,jlong ref,jlong destination){
  ((precice::ReceiverNativeDispatcher*)ref)->connect((precice::Receiver*)destination);
}

JNIEXPORT void JNICALL Java_precice_ReceiverNativeDispatcher_disconnect(JNIEnv *env, jobject obj,jlong ref,jlong destination){
  ((precice::ReceiverNativeDispatcher*)ref)->disconnect((precice::Receiver*)destination);
}


precice::ReceiverNativeDispatcher::ReceiverNativeDispatcher(){

}

precice::ReceiverNativeDispatcher::~ReceiverNativeDispatcher(){

}

void precice::ReceiverNativeDispatcher::connect(precice::Receiver* destination){
  if(std::find(_destinations.begin(), _destinations.end(), destination)==_destinations.end())
     _destinations.push_back(destination);
}

void precice::ReceiverNativeDispatcher::disconnect(precice::Receiver* destination){
  std::vector<precice::Receiver*>::iterator iter=std::find(_destinations.begin(), _destinations.end(), destination);
  if(iter!=_destinations.end())
     _destinations.erase(iter);
}

bool precice::ReceiverNativeDispatcher::isConnected() const{
  return !_destinations.empty();
}


void precice::ReceiverNativeDispatcher::receive(const double data,const int index,const int rank,int& tag){
    for(unsigned int i=0;i<_destinations.size();i++)
          _destinations[i]->receive(data,index,rank,tag);
}

void precice::ReceiverNativeDispatcher::receiveParallel(const double data,const int index,const int rank,int& tag){
    for(unsigned int i=0;i<_destinations.size();i++)
          _destinations[i]->receiveParallel(data,index,rank,tag);
}

