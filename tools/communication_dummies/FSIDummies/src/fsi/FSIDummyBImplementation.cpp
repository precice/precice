#include "fsi/FSIDummyBImplementation.h"

fsi::FSIDummyBImplementation::FSIDummyBImplementation(){

}

fsi::FSIDummyBImplementation::~FSIDummyBImplementation(){

}

extern "C" void
#ifdef _WIN32
MAIN_LOOP(bool joinable);
#else
main_loop_(bool joinable);
#endif

int main(){
#ifdef _WIN32
  MAIN_LOOP(false);
#else
  main_loop_(false);
#endif
}

void fsi::FSIDummyBImplementation::transferCoordinates(const double* coord, const int coord_len){
     std::cout<<"start receiving coord"<<std::endl;
     for(int i=0;i<coord_len;i++){
    	 std::cout<<"coord["<<i<<"]:"<<coord[i]<<std::endl;

     }
     std::cout<<"end receiving coord"<<std::endl;

}
void fsi::FSIDummyBImplementation::transferData(const double* data, const int data_len){
     // @todo Insert your code here
}
