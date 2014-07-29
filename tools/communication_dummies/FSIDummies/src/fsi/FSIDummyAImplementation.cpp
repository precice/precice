#include "fsi/FSIDummyAImplementation.h"
#include <iostream>

fsi::FSIDummyAImplementation::FSIDummyAImplementation(){
	std::cout<<"initializing FSIA"<<std::endl;
}

fsi::FSIDummyAImplementation::~FSIDummyAImplementation(){

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

void fsi::FSIDummyAImplementation::transferAllData(){
	double coordinates [2];
	coordinates[0]=0.0;
	coordinates[1]=0.0;
	_b->transferCoordinates(coordinates,2);
}

void fsi::FSIDummyAImplementation::test(){
	transferAllData();
}
