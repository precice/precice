#include "fsi/FSITestNativeDispatcher.h"
#include <algorithm>

JNIEXPORT void JNICALL Java_fsi_FSITestNativeDispatcher_createInstance(JNIEnv *env, jobject obj){
  JavaVM* jvm;
  env->GetJavaVM(&jvm);
  
  fsi::FSITestNativeDispatcher *ref=new fsi::FSITestNativeDispatcher();
  
  jfieldID id =env->GetFieldID(env->GetObjectClass(obj), "_ref", "J");
  env->SetLongField(obj, id, (jlong)ref);
    
}

JNIEXPORT void JNICALL Java_fsi_FSITestNativeDispatcher_destroyInstance(JNIEnv *env, jobject obj,jlong ref){
  delete ((fsi::FSITestNativeDispatcher*)ref);
}

JNIEXPORT void JNICALL Java_fsi_FSITestNativeDispatcher_connect(JNIEnv *env, jobject obj,jlong ref,jlong destination){
  ((fsi::FSITestNativeDispatcher*)ref)->connect((fsi::FSITest*)destination);
}

JNIEXPORT void JNICALL Java_fsi_FSITestNativeDispatcher_disconnect(JNIEnv *env, jobject obj,jlong ref,jlong destination){
  ((fsi::FSITestNativeDispatcher*)ref)->disconnect((fsi::FSITest*)destination);
}


fsi::FSITestNativeDispatcher::FSITestNativeDispatcher(){

}

fsi::FSITestNativeDispatcher::~FSITestNativeDispatcher(){

}

void fsi::FSITestNativeDispatcher::connect(fsi::FSITest* destination){
  if(std::find(_destinations.begin(), _destinations.end(), destination)==_destinations.end())
     _destinations.push_back(destination);
}

void fsi::FSITestNativeDispatcher::disconnect(fsi::FSITest* destination){
  std::vector<fsi::FSITest*>::iterator iter=std::find(_destinations.begin(), _destinations.end(), destination);
  if(iter!=_destinations.end())
     _destinations.erase(iter);
}

bool fsi::FSITestNativeDispatcher::isConnected() const{
  return !_destinations.empty();
}


void fsi::FSITestNativeDispatcher::test(){
    for(unsigned int i=0;i<_destinations.size();i++)
          _destinations[i]->test();
}

void fsi::FSITestNativeDispatcher::testParallel(){
    for(unsigned int i=0;i<_destinations.size();i++)
          _destinations[i]->testParallel();
}

