#include "precice/InitializerNativeDispatcher.h"
#include <algorithm>

JNIEXPORT void JNICALL Java_precice_InitializerNativeDispatcher_createInstance(JNIEnv *env, jobject obj){
  JavaVM* jvm;
  env->GetJavaVM(&jvm);
  
  precice::InitializerNativeDispatcher *ref=new precice::InitializerNativeDispatcher();
  
  jfieldID id =env->GetFieldID(env->GetObjectClass(obj), "_ref", "J");
  env->SetLongField(obj, id, (jlong)ref);
    
}

JNIEXPORT void JNICALL Java_precice_InitializerNativeDispatcher_destroyInstance(JNIEnv *env, jobject obj,jlong ref){
  delete ((precice::InitializerNativeDispatcher*)ref);
}

JNIEXPORT void JNICALL Java_precice_InitializerNativeDispatcher_connect(JNIEnv *env, jobject obj,jlong ref,jlong destination){
  ((precice::InitializerNativeDispatcher*)ref)->connect((precice::Initializer*)destination);
}

JNIEXPORT void JNICALL Java_precice_InitializerNativeDispatcher_disconnect(JNIEnv *env, jobject obj,jlong ref,jlong destination){
  ((precice::InitializerNativeDispatcher*)ref)->disconnect((precice::Initializer*)destination);
}


precice::InitializerNativeDispatcher::InitializerNativeDispatcher(){

}

precice::InitializerNativeDispatcher::~InitializerNativeDispatcher(){

}

void precice::InitializerNativeDispatcher::connect(precice::Initializer* destination){
  if(std::find(_destinations.begin(), _destinations.end(), destination)==_destinations.end())
     _destinations.push_back(destination);
}

void precice::InitializerNativeDispatcher::disconnect(precice::Initializer* destination){
  std::vector<precice::Initializer*>::iterator iter=std::find(_destinations.begin(), _destinations.end(), destination);
  if(iter!=_destinations.end())
     _destinations.erase(iter);
}

bool precice::InitializerNativeDispatcher::isConnected() const{
  return !_destinations.empty();
}


void precice::InitializerNativeDispatcher::acknowledge(const int identifier,int& tag){
    for(unsigned int i=0;i<_destinations.size();i++)
          _destinations[i]->acknowledge(identifier,tag);
}

void precice::InitializerNativeDispatcher::acknowledgeParallel(const int identifier,int& tag){
    for(unsigned int i=0;i<_destinations.size();i++)
          _destinations[i]->acknowledgeParallel(identifier,tag);
}
void precice::InitializerNativeDispatcher::initialize(const std::string* addresses, const int addresses_len,const int* vertexes, const int vertexes_len){
    for(unsigned int i=0;i<_destinations.size();i++)
          _destinations[i]->initialize(addresses,addresses_len,vertexes,vertexes_len);
}

void precice::InitializerNativeDispatcher::initializeParallel(const std::string* addresses, const int addresses_len,const int* vertexes, const int vertexes_len){
    for(unsigned int i=0;i<_destinations.size();i++)
          _destinations[i]->initializeParallel(addresses,addresses_len,vertexes,vertexes_len);
}

