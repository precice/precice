#include "fsi/FSIDummyBAbstractImplementation.h"

fsi::FSIDummyBAbstractImplementation::FSIDummyBAbstractImplementation(){
     
}

fsi::FSIDummyBAbstractImplementation::~FSIDummyBAbstractImplementation(){

}

void fsi::FSIDummyBAbstractImplementation::transferDataParallel(const double* data, const int data_len){
     // @todo Insert your code here
}
void fsi::FSIDummyBAbstractImplementation::transferCoordinatesParallel(const int* coordId, const int coordId_len,const int* offsets, const int offsets_len,const std::string* hosts, const int hosts_len){
     // @todo Insert your code here
}
