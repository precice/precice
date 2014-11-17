#include "precice/ReceiverNative2JavaPlainPort.h"

JNIEXPORT void JNICALL Java_precice_ReceiverNative2JavaPlainPort_createInstance(JNIEnv *env, jobject obj){
  JavaVM* jvm;
  env->GetJavaVM(&jvm);
  jobject self=env->NewGlobalRef(obj);
  
  precice::ReceiverNative2JavaPlainPort *ref=new precice::ReceiverNative2JavaPlainPort(jvm,self);
  
  jfieldID id =env->GetFieldID(env->GetObjectClass(obj), "_ref", "J");
  env->SetLongField(obj, id, (jlong)ref);
  
}

JNIEXPORT void JNICALL Java_precice_ReceiverNative2JavaPlainPort_destroyInstance(JNIEnv *env, jobject obj,jlong ref){
  delete ((precice::ReceiverNative2JavaPlainPort*)ref);
  
}

precice::ReceiverNative2JavaPlainPort::ReceiverNative2JavaPlainPort(JavaVM* jvm,jobject obj):
     _jvm(jvm),
     _obj(obj){

}

precice::ReceiverNative2JavaPlainPort::~ReceiverNative2JavaPlainPort(){
  JNIEnv* env;
  _jvm->GetEnv((void**)&env,JNI_VERSION_1_6);
  env->DeleteGlobalRef(_obj);
}

void precice::ReceiverNative2JavaPlainPort::receive(const double data,const int index,const int rank,int& tag){
  JNIEnv* env;
  int status=_jvm->GetEnv((void**)&env,JNI_VERSION_1_6);
  jclass cls=env->GetObjectClass(_obj);
  jmethodID mid = env->GetMethodID(cls,"receive","(DII[I)V");
  
  if(mid){
     jdouble data_jni=data;
jint index_jni=index;
jint rank_jni=rank;
jintArray tag_jni=env->NewIntArray(1);
env->SetIntArrayRegion(tag_jni,0,1,(jint*)&tag);

     env->CallVoidMethod(_obj,mid,data_jni,index_jni,rank_jni,tag_jni);
     env->GetIntArrayRegion(tag_jni,0,1,(jint*)&tag);

  }
}
void precice::ReceiverNative2JavaPlainPort::receiveParallel(const double data,const int index,const int rank,int& tag){
  JNIEnv* env;
  int status=_jvm->GetEnv((void**)&env,JNI_VERSION_1_6);
  jclass cls=env->GetObjectClass(_obj);
  jmethodID mid = env->GetMethodID(cls,"receiveParallel","(DII[I)V");
  
  if(mid){
     jdouble data_jni=data;
jint index_jni=index;
jint rank_jni=rank;
jintArray tag_jni=env->NewIntArray(1);
env->SetIntArrayRegion(tag_jni,0,1,(jint*)&tag);

     env->CallVoidMethod(_obj,mid,data_jni,index_jni,rank_jni,tag_jni);
     env->GetIntArrayRegion(tag_jni,0,1,(jint*)&tag);

  }
}
