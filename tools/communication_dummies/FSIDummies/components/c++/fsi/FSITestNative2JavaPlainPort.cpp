#include "fsi/FSITestNative2JavaPlainPort.h"

JNIEXPORT void JNICALL Java_fsi_FSITestNative2JavaPlainPort_createInstance(JNIEnv *env, jobject obj){
  JavaVM* jvm;
  env->GetJavaVM(&jvm);
  jobject self=env->NewGlobalRef(obj);
  
  fsi::FSITestNative2JavaPlainPort *ref=new fsi::FSITestNative2JavaPlainPort(jvm,self);
  
  jfieldID id =env->GetFieldID(env->GetObjectClass(obj), "_ref", "J");
  env->SetLongField(obj, id, (jlong)ref);
  
}

JNIEXPORT void JNICALL Java_fsi_FSITestNative2JavaPlainPort_destroyInstance(JNIEnv *env, jobject obj,jlong ref){
  delete ((fsi::FSITestNative2JavaPlainPort*)ref);
  
}

fsi::FSITestNative2JavaPlainPort::FSITestNative2JavaPlainPort(JavaVM* jvm,jobject obj):
     _jvm(jvm),
     _obj(obj){

}

fsi::FSITestNative2JavaPlainPort::~FSITestNative2JavaPlainPort(){
  JNIEnv* env;
  _jvm->GetEnv((void**)&env,JNI_VERSION_1_6);
  env->DeleteGlobalRef(_obj);
}

void fsi::FSITestNative2JavaPlainPort::test(){
  JNIEnv* env;
  int status=_jvm->GetEnv((void**)&env,JNI_VERSION_1_6);
  jclass cls=env->GetObjectClass(_obj);
  jmethodID mid = env->GetMethodID(cls,"test","()V");
  
  if(mid){
     
     env->CallVoidMethod(_obj,mid);
     
  }
}
void fsi::FSITestNative2JavaPlainPort::testParallel(){
  JNIEnv* env;
  int status=_jvm->GetEnv((void**)&env,JNI_VERSION_1_6);
  jclass cls=env->GetObjectClass(_obj);
  jmethodID mid = env->GetMethodID(cls,"testParallel","()V");
  
  if(mid){
     
     env->CallVoidMethod(_obj,mid);
     
  }
}
