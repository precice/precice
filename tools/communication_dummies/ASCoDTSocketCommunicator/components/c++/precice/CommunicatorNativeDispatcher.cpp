#include "precice/CommunicatorNativeDispatcher.h"
#include <algorithm>

JNIEXPORT void JNICALL Java_precice_CommunicatorNativeDispatcher_createInstance(JNIEnv *env, jobject obj){
  JavaVM* jvm;
  env->GetJavaVM(&jvm);
  
  precice::CommunicatorNativeDispatcher *ref=new precice::CommunicatorNativeDispatcher();
  
  jfieldID id =env->GetFieldID(env->GetObjectClass(obj), "_ref", "J");
  env->SetLongField(obj, id, (jlong)ref);
    
}

JNIEXPORT void JNICALL Java_precice_CommunicatorNativeDispatcher_destroyInstance(JNIEnv *env, jobject obj,jlong ref){
  delete ((precice::CommunicatorNativeDispatcher*)ref);
}

JNIEXPORT void JNICALL Java_precice_CommunicatorNativeDispatcher_connect(JNIEnv *env, jobject obj,jlong ref,jlong destination){
  ((precice::CommunicatorNativeDispatcher*)ref)->connect((precice::Communicator*)destination);
}

JNIEXPORT void JNICALL Java_precice_CommunicatorNativeDispatcher_disconnect(JNIEnv *env, jobject obj,jlong ref,jlong destination){
  ((precice::CommunicatorNativeDispatcher*)ref)->disconnect((precice::Communicator*)destination);
}


precice::CommunicatorNativeDispatcher::CommunicatorNativeDispatcher(){

}

precice::CommunicatorNativeDispatcher::~CommunicatorNativeDispatcher(){

}

void precice::CommunicatorNativeDispatcher::connect(precice::Communicator* destination){
  if(std::find(_destinations.begin(), _destinations.end(), destination)==_destinations.end())
     _destinations.push_back(destination);
}

void precice::CommunicatorNativeDispatcher::disconnect(precice::Communicator* destination){
  std::vector<precice::Communicator*>::iterator iter=std::find(_destinations.begin(), _destinations.end(), destination);
  if(iter!=_destinations.end())
     _destinations.erase(iter);
}

bool precice::CommunicatorNativeDispatcher::isConnected() const{
  return !_destinations.empty();
}


void precice::CommunicatorNativeDispatcher::setData(const double* data, const int data_len){
    for(unsigned int i=0;i<_destinations.size();i++)
          _destinations[i]->setData(data,data_len);
}

void precice::CommunicatorNativeDispatcher::setDataParallel(const double* data, const int data_len){
    for(unsigned int i=0;i<_destinations.size();i++)
          _destinations[i]->setDataParallel(data,data_len);
}

