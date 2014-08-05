#include "mpi.h"
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>
#include "math.h"

#include "fsi/FSIDummyBImplementation.h"
/**
 * @brief For printing to the command line.
 *
 * @param message  An input stream such as: "Blabla = " << variable << ", ..."
 */
#define PRINT(message) \
  { \
    std::ostringstream conv; \
    conv << message; \
    std::cout << conv.str() << std::endl; \
  }
extern "C" void main_loop_(bool);
int main(int argc, char** argv)
{
  std::cout << "Running communication proxy" << std::endl;
  int provided;
  MPI_Init_thread(&argc,&argv,MPI_THREAD_MULTIPLE ,&provided);

  if (argc == 1 || argc != 3){
    PRINT("Usage: ./proxy proX proY");
    return 1;
  }
  std::cout<<"argv2:"<<argv[2]<<std::endl;
  int proX = atoi(argv[1]);
  int proY = atoi(argv[2]);

  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (proX * proY != size){
    PRINT("proX * proY should be size");
    return 1;
  }

  int rankX = rank % proX +1;
  int rankY = rank / proX + 1;

  int N = 10;

  int chunkX, chunkY;

  if(rankX < proX){
    chunkX  = N / proX;
  }
  else{
    chunkX = N / proX + N % proX;
  }
  if(rankY < proY){
    chunkY  = N / proY;
  }
  else{
    chunkY = N / proY + N % proY;
  }

  int pointSize = chunkX*chunkY;
  double* coordX = new double[pointSize];
  double* coordY = new double[pointSize];
  double* data = new double[pointSize];
  int* ids = new int[pointSize];

  double offsetX = (rankX-1)* (N/proX);
  double offsetY = (rankY-1)* (N/proY);

 

  double dx = 1.0 / N;
  double dy = 1.0 / N;

  for(int i=1;i<=chunkX;i++){
    for(int j=1;j<=chunkY;j++){
      coordX[(j-1)*chunkX +i-1] = (i-1)*dx + offsetX*dx;
      coordY[(j-1)*chunkX +i-1] = (j-1)*dy + offsetY*dy;
    }
  }

  for(int k=0;k<pointSize;k++){
    ids[k] = round((coordX[k]/dx)+(coordY[k]/dy)*N);
  }
 


  PRINT("Second Step:B");
  //B sends data

  for(int k=0;k<pointSize;k++){
    data[k]=coordX[k]*coordY[k]*coordY[k]+1.0;
  }

  main_loop_(false);
  fsi::FSIDummyBImplementation::singleton->setCoordinates(&ids[0],pointSize);
  fsi::FSIDummyBImplementation::singleton->setData(data);
  int a = 20;
  while(true){
	a=a+30;

  }
  MPI_Finalize();
  std::cout << "Stop communication proxy" << std::endl;
  return 0;
}

