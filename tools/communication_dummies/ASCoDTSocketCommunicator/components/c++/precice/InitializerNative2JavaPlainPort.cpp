#include "precice/InitializerNative2JavaPlainPort.h"

JNIEXPORT void JNICALL Java_precice_InitializerNative2JavaPlainPort_createInstance(JNIEnv *env, jobject obj){
  JavaVM* jvm;
  env->GetJavaVM(&jvm);
  jobject self=env->NewGlobalRef(obj);
  
  precice::InitializerNative2JavaPlainPort *ref=new precice::InitializerNative2JavaPlainPort(jvm,self);
  
  jfieldID id =env->GetFieldID(env->GetObjectClass(obj), "_ref", "J");
  env->SetLongField(obj, id, (jlong)ref);
  
}

JNIEXPORT void JNICALL Java_precice_InitializerNative2JavaPlainPort_destroyInstance(JNIEnv *env, jobject obj,jlong ref){
  delete ((precice::InitializerNative2JavaPlainPort*)ref);
  
}

precice::InitializerNative2JavaPlainPort::InitializerNative2JavaPlainPort(JavaVM* jvm,jobject obj):
     _jvm(jvm),
     _obj(obj){

}

precice::InitializerNative2JavaPlainPort::~InitializerNative2JavaPlainPort(){
  JNIEnv* env;
  _jvm->GetEnv((void**)&env,JNI_VERSION_1_6);
  env->DeleteGlobalRef(_obj);
}

void precice::InitializerNative2JavaPlainPort::initializeAddresses(const std::string* addresses, const int addresses_len){
  JNIEnv* env;
  int status=_jvm->GetEnv((void**)&env,JNI_VERSION_1_6);
  jclass cls=env->GetObjectClass(_obj);
  jmethodID mid = env->GetMethodID(cls,"initializeAddresses","([Ljava/lang/String;)V");
  
  if(mid){
     jobjectArray addresses_jni=env->NewObjectArray(addresses_len,env->FindClass("Ljava/lang/String;"),0);
for(int i=0;i<addresses_len;i++)
	env->SetObjectArrayElement(addresses_jni,i, env->NewStringUTF(addresses[i].c_str()));

     env->CallVoidMethod(_obj,mid,addresses_jni);
     
  }
}
void precice::InitializerNative2JavaPlainPort::initializeAddressesParallel(const std::string* addresses, const int addresses_len){
  JNIEnv* env;
  int status=_jvm->GetEnv((void**)&env,JNI_VERSION_1_6);
  jclass cls=env->GetObjectClass(_obj);
  jmethodID mid = env->GetMethodID(cls,"initializeAddressesParallel","([Ljava/lang/String;)V");
  
  if(mid){
     jobjectArray addresses_jni=env->NewObjectArray(addresses_len,env->FindClass("Ljava/lang/String;"),0);
for(int i=0;i<addresses_len;i++)
	env->SetObjectArrayElement(addresses_jni,i, env->NewStringUTF(addresses[i].c_str()));

     env->CallVoidMethod(_obj,mid,addresses_jni);
     
  }
}
void precice::InitializerNative2JavaPlainPort::initializeVertexes(const int* vertexes, const int vertexes_len){
  JNIEnv* env;
  int status=_jvm->GetEnv((void**)&env,JNI_VERSION_1_6);
  jclass cls=env->GetObjectClass(_obj);
  jmethodID mid = env->GetMethodID(cls,"initializeVertexes","([I)V");
  
  if(mid){
     jintArray vertexes_jni=env->NewIntArray(vertexes_len);
env->SetIntArrayRegion(vertexes_jni,0,vertexes_len,(jint*)vertexes);

     env->CallVoidMethod(_obj,mid,vertexes_jni);
     
  }
}
void precice::InitializerNative2JavaPlainPort::initializeVertexesParallel(const int* vertexes, const int vertexes_len){
  JNIEnv* env;
  int status=_jvm->GetEnv((void**)&env,JNI_VERSION_1_6);
  jclass cls=env->GetObjectClass(_obj);
  jmethodID mid = env->GetMethodID(cls,"initializeVertexesParallel","([I)V");
  
  if(mid){
     jintArray vertexes_jni=env->NewIntArray(vertexes_len);
env->SetIntArrayRegion(vertexes_jni,0,vertexes_len,(jint*)vertexes);

     env->CallVoidMethod(_obj,mid,vertexes_jni);
     
  }
}
