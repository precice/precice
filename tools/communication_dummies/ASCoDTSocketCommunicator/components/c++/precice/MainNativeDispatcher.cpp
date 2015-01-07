#include "precice/MainNativeDispatcher.h"
#include <algorithm>

JNIEXPORT void JNICALL Java_precice_MainNativeDispatcher_createInstance(JNIEnv *env, jobject obj){
  JavaVM* jvm;
  env->GetJavaVM(&jvm);
  
  precice::MainNativeDispatcher *ref=new precice::MainNativeDispatcher();
  
  jfieldID id =env->GetFieldID(env->GetObjectClass(obj), "_ref", "J");
  env->SetLongField(obj, id, (jlong)ref);
    
}

JNIEXPORT void JNICALL Java_precice_MainNativeDispatcher_destroyInstance(JNIEnv *env, jobject obj,jlong ref){
  delete ((precice::MainNativeDispatcher*)ref);
}

JNIEXPORT void JNICALL Java_precice_MainNativeDispatcher_connect(JNIEnv *env, jobject obj,jlong ref,jlong destination){
  ((precice::MainNativeDispatcher*)ref)->connect((precice::Main*)destination);
}

JNIEXPORT void JNICALL Java_precice_MainNativeDispatcher_disconnect(JNIEnv *env, jobject obj,jlong ref,jlong destination){
  ((precice::MainNativeDispatcher*)ref)->disconnect((precice::Main*)destination);
}


precice::MainNativeDispatcher::MainNativeDispatcher(){

}

precice::MainNativeDispatcher::~MainNativeDispatcher(){

}

void precice::MainNativeDispatcher::connect(precice::Main* destination){
  if(std::find(_destinations.begin(), _destinations.end(), destination)==_destinations.end())
     _destinations.push_back(destination);
}

void precice::MainNativeDispatcher::disconnect(precice::Main* destination){
  std::vector<precice::Main*>::iterator iter=std::find(_destinations.begin(), _destinations.end(), destination);
  if(iter!=_destinations.end())
     _destinations.erase(iter);
}

bool precice::MainNativeDispatcher::isConnected() const{
  return !_destinations.empty();
}


void precice::MainNativeDispatcher::main(){
    for(unsigned int i=0;i<_destinations.size();i++)
          _destinations[i]->main();
}

void precice::MainNativeDispatcher::mainParallel(){
    for(unsigned int i=0;i<_destinations.size();i++)
          _destinations[i]->mainParallel();
}

