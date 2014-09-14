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

void precice::CommunicatorNative2JavaPlainPort::setData(const double* data, const int data_len){
  JNIEnv* env;
  int status=_jvm->GetEnv((void**)&env,JNI_VERSION_1_6);
  jclass cls=env->GetObjectClass(_obj);
  jmethodID mid = env->GetMethodID(cls,"setData","([D)V");
  
  if(mid){
     jdoubleArray data_jni=env->NewDoubleArray(data_len);
env->SetDoubleArrayRegion(data_jni,0,data_len,(jdouble*)data);

     env->CallVoidMethod(_obj,mid,data_jni);
     
  }
}
void precice::CommunicatorNative2JavaPlainPort::setDataParallel(const double* data, const int data_len){
  JNIEnv* env;
  int status=_jvm->GetEnv((void**)&env,JNI_VERSION_1_6);
  jclass cls=env->GetObjectClass(_obj);
  jmethodID mid = env->GetMethodID(cls,"setDataParallel","([D)V");
  
  if(mid){
     jdoubleArray data_jni=env->NewDoubleArray(data_len);
env->SetDoubleArrayRegion(data_jni,0,data_len,(jdouble*)data);

     env->CallVoidMethod(_obj,mid,data_jni);
     
  }
}
