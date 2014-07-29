#include "fsi/FSICommNative2JavaPlainPort.h"

JNIEXPORT void JNICALL Java_fsi_FSICommNative2JavaPlainPort_createInstance(JNIEnv *env, jobject obj){
  JavaVM* jvm;
  env->GetJavaVM(&jvm);
  jobject self=env->NewGlobalRef(obj);
  
  fsi::FSICommNative2JavaPlainPort *ref=new fsi::FSICommNative2JavaPlainPort(jvm,self);
  
  jfieldID id =env->GetFieldID(env->GetObjectClass(obj), "_ref", "J");
  env->SetLongField(obj, id, (jlong)ref);
  
}

JNIEXPORT void JNICALL Java_fsi_FSICommNative2JavaPlainPort_destroyInstance(JNIEnv *env, jobject obj,jlong ref){
  delete ((fsi::FSICommNative2JavaPlainPort*)ref);
  
}

fsi::FSICommNative2JavaPlainPort::FSICommNative2JavaPlainPort(JavaVM* jvm,jobject obj):
     _jvm(jvm),
     _obj(obj){

}

fsi::FSICommNative2JavaPlainPort::~FSICommNative2JavaPlainPort(){
  JNIEnv* env;
  _jvm->GetEnv((void**)&env,JNI_VERSION_1_6);
  env->DeleteGlobalRef(_obj);
}

void fsi::FSICommNative2JavaPlainPort::transferCoordinates(const double* coord, const int coord_len){
  JNIEnv* env;
  int status=_jvm->GetEnv((void**)&env,JNI_VERSION_1_6);
  jclass cls=env->GetObjectClass(_obj);
  jmethodID mid = env->GetMethodID(cls,"transferCoordinates","([D)V");
  
  if(mid){
     jdoubleArray coord_jni=env->NewDoubleArray(coord_len);
env->SetDoubleArrayRegion(coord_jni,0,coord_len,(jdouble*)coord);

     env->CallVoidMethod(_obj,mid,coord_jni);
     
  }
}
void fsi::FSICommNative2JavaPlainPort::transferCoordinatesParallel(const double* coord, const int coord_len){
  JNIEnv* env;
  int status=_jvm->GetEnv((void**)&env,JNI_VERSION_1_6);
  jclass cls=env->GetObjectClass(_obj);
  jmethodID mid = env->GetMethodID(cls,"transferCoordinatesParallel","([D)V");
  
  if(mid){
     jdoubleArray coord_jni=env->NewDoubleArray(coord_len);
env->SetDoubleArrayRegion(coord_jni,0,coord_len,(jdouble*)coord);

     env->CallVoidMethod(_obj,mid,coord_jni);
     
  }
}
void fsi::FSICommNative2JavaPlainPort::transferData(const double* data, const int data_len){
  JNIEnv* env;
  int status=_jvm->GetEnv((void**)&env,JNI_VERSION_1_6);
  jclass cls=env->GetObjectClass(_obj);
  jmethodID mid = env->GetMethodID(cls,"transferData","([D)V");
  
  if(mid){
     jdoubleArray data_jni=env->NewDoubleArray(data_len);
env->SetDoubleArrayRegion(data_jni,0,data_len,(jdouble*)data);

     env->CallVoidMethod(_obj,mid,data_jni);
     
  }
}
void fsi::FSICommNative2JavaPlainPort::transferDataParallel(const double* data, const int data_len){
  JNIEnv* env;
  int status=_jvm->GetEnv((void**)&env,JNI_VERSION_1_6);
  jclass cls=env->GetObjectClass(_obj);
  jmethodID mid = env->GetMethodID(cls,"transferDataParallel","([D)V");
  
  if(mid){
     jdoubleArray data_jni=env->NewDoubleArray(data_len);
env->SetDoubleArrayRegion(data_jni,0,data_len,(jdouble*)data);

     env->CallVoidMethod(_obj,mid,data_jni);
     
  }
}
