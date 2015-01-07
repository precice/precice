#include "fsi/FSICommNative2JavaPlainPort.h"

#ifdef JAVA

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

void fsi::FSICommNative2JavaPlainPort::transferCoordinates(const int* coordId, const int coordId_len,const int* offsets, const int offsets_len,const std::string* hosts, const int hosts_len){
  JNIEnv* env;
  int status=_jvm->GetEnv((void**)&env,JNI_VERSION_1_6);
  jclass cls=env->GetObjectClass(_obj);
  jmethodID mid = env->GetMethodID(cls,"transferCoordinates","([I[I[Ljava/lang/String;)V");
  
  if(mid){
     jintArray coordId_jni=env->NewIntArray(coordId_len);
env->SetIntArrayRegion(coordId_jni,0,coordId_len,(jint*)coordId);
jintArray offsets_jni=env->NewIntArray(offsets_len);
env->SetIntArrayRegion(offsets_jni,0,offsets_len,(jint*)offsets);
jobjectArray hosts_jni=env->NewObjectArray(hosts_len,env->FindClass("Ljava/lang/String;"),0);
for(int i=0;i<hosts_len;i++)
	env->SetObjectArrayElement(hosts_jni,i, env->NewStringUTF(hosts[i].c_str()));

     env->CallVoidMethod(_obj,mid,coordId_jni,offsets_jni,hosts_jni);
     
  }
}
void fsi::FSICommNative2JavaPlainPort::transferCoordinatesParallel(const int* coordId, const int coordId_len,const int* offsets, const int offsets_len,const std::string* hosts, const int hosts_len){
  JNIEnv* env;
  int status=_jvm->GetEnv((void**)&env,JNI_VERSION_1_6);
  jclass cls=env->GetObjectClass(_obj);
  jmethodID mid = env->GetMethodID(cls,"transferCoordinatesParallel","([I[I[Ljava/lang/String;)V");
  
  if(mid){
     jintArray coordId_jni=env->NewIntArray(coordId_len);
env->SetIntArrayRegion(coordId_jni,0,coordId_len,(jint*)coordId);
jintArray offsets_jni=env->NewIntArray(offsets_len);
env->SetIntArrayRegion(offsets_jni,0,offsets_len,(jint*)offsets);
jobjectArray hosts_jni=env->NewObjectArray(hosts_len,env->FindClass("Ljava/lang/String;"),0);
for(int i=0;i<hosts_len;i++)
	env->SetObjectArrayElement(hosts_jni,i, env->NewStringUTF(hosts[i].c_str()));

     env->CallVoidMethod(_obj,mid,coordId_jni,offsets_jni,hosts_jni);
     
  }
}
void fsi::FSICommNative2JavaPlainPort::startDataTransfer(){
  JNIEnv* env;
  int status=_jvm->GetEnv((void**)&env,JNI_VERSION_1_6);
  jclass cls=env->GetObjectClass(_obj);
  jmethodID mid = env->GetMethodID(cls,"startDataTransfer","()V");
  
  if(mid){
     
     env->CallVoidMethod(_obj,mid);
     
  }
}
void fsi::FSICommNative2JavaPlainPort::startDataTransferParallel(){
  JNIEnv* env;
  int status=_jvm->GetEnv((void**)&env,JNI_VERSION_1_6);
  jclass cls=env->GetObjectClass(_obj);
  jmethodID mid = env->GetMethodID(cls,"startDataTransferParallel","()V");
  
  if(mid){
     
     env->CallVoidMethod(_obj,mid);
     
  }
}
void fsi::FSICommNative2JavaPlainPort::endDataTransfer(int& ack){
  JNIEnv* env;
  int status=_jvm->GetEnv((void**)&env,JNI_VERSION_1_6);
  jclass cls=env->GetObjectClass(_obj);
  jmethodID mid = env->GetMethodID(cls,"endDataTransfer","([I)V");
  
  if(mid){
     jintArray ack_jni=env->NewIntArray(1);
env->SetIntArrayRegion(ack_jni,0,1,(jint*)&ack);

     env->CallVoidMethod(_obj,mid,ack_jni);
     env->GetIntArrayRegion(ack_jni,0,1,(jint*)&ack);

  }
}
void fsi::FSICommNative2JavaPlainPort::endDataTransferParallel(int& ack){
  JNIEnv* env;
  int status=_jvm->GetEnv((void**)&env,JNI_VERSION_1_6);
  jclass cls=env->GetObjectClass(_obj);
  jmethodID mid = env->GetMethodID(cls,"endDataTransferParallel","([I)V");
  
  if(mid){
     jintArray ack_jni=env->NewIntArray(1);
env->SetIntArrayRegion(ack_jni,0,1,(jint*)&ack);

     env->CallVoidMethod(_obj,mid,ack_jni);
     env->GetIntArrayRegion(ack_jni,0,1,(jint*)&ack);

  }
}
#endif
