#include "precice/BAbstractImplementation.h"

precice::BAbstractImplementation::BAbstractImplementation(){
     _a = 0;

}

precice::BAbstractImplementation::~BAbstractImplementation(){

}

/**
 * @see Case class 
 */
 void precice::BAbstractImplementation::connecta(precice::InitializerNativeDispatcher* port){
   _a = port; 
 }
 
 void precice::BAbstractImplementation::disconnecta(){
    //delete _a;
    _a=NULL;
 }
 
void precice::BAbstractImplementation::setDataParallel(const double data,const int index,const int rank,int& tag){
     // @todo Insert your code here
}
void precice::BAbstractImplementation::initializeParallel(const std::string* addresses, const int addresses_len,const int* vertexes, const int vertexes_len){
     // @todo Insert your code here
}
void precice::BAbstractImplementation::acknowledgeParallel(const int identifier,int& tag){
     // @todo Insert your code here
}
void precice::BAbstractImplementation::mainParallel(){
     // @todo Insert your code here
}
