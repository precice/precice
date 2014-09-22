#include "precice/CommunicatorNative2JavaPlainPort.h"

JNIEXPORT void JNICALL Java_precice_CommunicatorNative2JavaPlainPort_createInstance(JNIEnv *env, jobject obj){
  JavaVM* jvm;
  env->GetJavaVM(&jvm);
  jobject self=env->NewGlobalRef(obj);
  
  precice::CommunicatorNative2JavaPlainPort *ref=new precice::CommunicatorNative2JavaPlainPort(jvm,self);
  
  jfieldID id =env->GetFieldID(env->GetObjectClass(obj), "_ref", "J");
  env->SetLongField(obj, id, (jlong)ref);
  
}

JNIEXPORT void JNICALL Java_precice_CommunicatorNative2JavaPlainPort_destroyInstance(JNIEnv *env, jobject obj,jlong ref){
  delete ((precice::CommunicatorNative2JavaPlainPort*)ref);
  
}

precice::CommunicatorNative2JavaPlainPort::CommunicatorNative2JavaPlainPort(JavaVM* jvm,jobject obj):
     _jvm(jvm),
     _obj(obj){

}

precice::CommunicatorNative2JavaPlainPort::~CommunicatorNative2JavaPlainPort(){
  JNIEnv* env;
  _jvm->GetEnv((void**)&env,JNI_VERSION_1_6);
  env->DeleteGlobalRef(_obj);
}

void precice::CommunicatorNative2JavaPlainPort::setData(const double data,const int index,const int rank,int& tag){
  JNIEnv* env;
  int status=_jvm->GetEnv((void**)&env,JNI_VERSION_1_6);
  jclass cls=env->GetObjectClass(_obj);
  jmethodID mid = env->GetMethodID(cls,"setData","(DII[I)V");
  
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
void precice::CommunicatorNative2JavaPlainPort::setDataParallel(const double data,const int index,const int rank,int& tag){
  JNIEnv* env;
  int status=_jvm->GetEnv((void**)&env,JNI_VERSION_1_6);
  jclass cls=env->GetObjectClass(_obj);
  jmethodID mid = env->GetMethodID(cls,"setDataParallel","(DII[I)V");
  
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
