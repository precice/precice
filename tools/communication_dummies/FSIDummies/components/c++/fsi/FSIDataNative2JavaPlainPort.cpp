#include "fsi/FSIDataNative2JavaPlainPort.h"

#ifdef JAVA

JNIEXPORT void JNICALL Java_fsi_FSIDataNative2JavaPlainPort_createInstance(JNIEnv *env, jobject obj){
  JavaVM* jvm;
  env->GetJavaVM(&jvm);
  jobject self=env->NewGlobalRef(obj);
  
  fsi::FSIDataNative2JavaPlainPort *ref=new fsi::FSIDataNative2JavaPlainPort(jvm,self);
  
  jfieldID id =env->GetFieldID(env->GetObjectClass(obj), "_ref", "J");
  env->SetLongField(obj, id, (jlong)ref);
  
}

JNIEXPORT void JNICALL Java_fsi_FSIDataNative2JavaPlainPort_destroyInstance(JNIEnv *env, jobject obj,jlong ref){
  delete ((fsi::FSIDataNative2JavaPlainPort*)ref);
  
}



fsi::FSIDataNative2JavaPlainPort::FSIDataNative2JavaPlainPort(JavaVM* jvm,jobject obj):
     _jvm(jvm),
     _obj(obj){

}

fsi::FSIDataNative2JavaPlainPort::~FSIDataNative2JavaPlainPort(){
  JNIEnv* env;
  _jvm->GetEnv((void**)&env,JNI_VERSION_1_6);
  env->DeleteGlobalRef(_obj);
}

void fsi::FSIDataNative2JavaPlainPort::transferData(const int* coordId, const int coordId_len,const double* data, const int data_len){
  JNIEnv* env;
  int status=_jvm->GetEnv((void**)&env,JNI_VERSION_1_6);
  jclass cls=env->GetObjectClass(_obj);
  jmethodID mid = env->GetMethodID(cls,"transferData","([I[D)V");
  
  if(mid){
     jintArray coordId_jni=env->NewIntArray(coordId_len);
env->SetIntArrayRegion(coordId_jni,0,coordId_len,(jint*)coordId);
jdoubleArray data_jni=env->NewDoubleArray(data_len);
env->SetDoubleArrayRegion(data_jni,0,data_len,(jdouble*)data);

     env->CallVoidMethod(_obj,mid,coordId_jni,data_jni);
     
  }
}
void fsi::FSIDataNative2JavaPlainPort::transferDataParallel(const int* coordId, const int coordId_len,const double* data, const int data_len){
  JNIEnv* env;
  int status=_jvm->GetEnv((void**)&env,JNI_VERSION_1_6);
  jclass cls=env->GetObjectClass(_obj);
  jmethodID mid = env->GetMethodID(cls,"transferDataParallel","([I[D)V");
  
  if(mid){
     jintArray coordId_jni=env->NewIntArray(coordId_len);
env->SetIntArrayRegion(coordId_jni,0,coordId_len,(jint*)coordId);
jdoubleArray data_jni=env->NewDoubleArray(data_len);
env->SetDoubleArrayRegion(data_jni,0,data_len,(jdouble*)data);

     env->CallVoidMethod(_obj,mid,coordId_jni,data_jni);
     
  }
}
void fsi::FSIDataNative2JavaPlainPort::dataAck(int& ack){
  JNIEnv* env;
  int status=_jvm->GetEnv((void**)&env,JNI_VERSION_1_6);
  jclass cls=env->GetObjectClass(_obj);
  jmethodID mid = env->GetMethodID(cls,"dataAck","([I)V");
  
  if(mid){
     jintArray ack_jni=env->NewIntArray(1);
env->SetIntArrayRegion(ack_jni,0,1,(jint*)&ack);

     env->CallVoidMethod(_obj,mid,ack_jni);
     env->GetIntArrayRegion(ack_jni,0,1,(jint*)&ack);

  }
}
void fsi::FSIDataNative2JavaPlainPort::dataAckParallel(int& ack){
  JNIEnv* env;
  int status=_jvm->GetEnv((void**)&env,JNI_VERSION_1_6);
  jclass cls=env->GetObjectClass(_obj);
  jmethodID mid = env->GetMethodID(cls,"dataAckParallel","([I)V");
  
  if(mid){
     jintArray ack_jni=env->NewIntArray(1);
env->SetIntArrayRegion(ack_jni,0,1,(jint*)&ack);

     env->CallVoidMethod(_obj,mid,ack_jni);
     env->GetIntArrayRegion(ack_jni,0,1,(jint*)&ack);

  }
}
#endif
