#include "mpi.h"
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>
#include "math.h"

#include "fsi/FSIDummyAImplementation.h"

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
extern "C" struct FSI_FSIDUMMYA_arg;
extern "C" FSI_FSIDUMMYA_arg daemon_args;
extern "C" void main_loop_(bool);
extern "C" void initialise_(FSI_FSIDUMMYA_arg& arg,bool joinable);

extern "C" void bind_component_(FSI_FSIDUMMYA_arg& arg,bool joinable);
extern "C" void socket_loop_(FSI_FSIDUMMYA_arg& arg,bool joinable);
extern "C" void destroy_(FSI_FSIDUMMYA_arg& arg,bool joinable);  
int main(int argc, char** argv)
{
  std::cout << "Running communication proxy" << std::endl;
  int provided;
  MPI_Init_thread(&argc,&argv,MPI_THREAD_MULTIPLE ,&provided);

  if (argc == 1 || argc != 3){
    PRINT("Usage: ./proxy proX proY");
    return 1;
  }

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

  std::cout << "Servus ..." << rank << " of " << size;
  std::cout << " (" << rankX << "," << rankY << ")";
  std::cout << " (" << chunkX << "," << chunkY << ")";
  std::cout << " (" << offsetX << "," << offsetY << ")" <<std::endl;

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
  for(int k=0;k<pointSize;k++){
    std::cout << ids[k] <<", " << std::endl;
  }


  PRINT("First Step:A");
  //A sets coordinates and receives data

  for(int k=0;k<pointSize;k++){
    data[k]=0.0;
  }
//
std::cout<<"rank:"<<rank<<"init cmp"<<std::endl;  
initialise_(daemon_args,false);
  fsi::FSIDummyAImplementation::singleton->setCoordIds(ids,pointSize);
  fsi::FSIDummyAImplementation::singleton->setData(data);
  fsi::FSIDummyAImplementation::singleton->gatherMids();
  fsi::FSIDummyAImplementation::singleton->gatherDomainDescriptions();
  MPI_Barrier(MPI_COMM_WORLD);
  std::cout<<"rank:"<<rank<<"bind cmp"<<std::endl;
  bind_component_(daemon_args,false);
  socket_loop_(daemon_args,false);
  
//
  fsi::FSIDummyAImplementation::singleton->transferGlobalIds();
  fsi::FSIDummyAImplementation::singleton->receiveAllData();

  //check if received data is correct
  std::cout << "Performaning validation" << std::endl;
  bool hasFailed=false;
  for(int k=0;k<pointSize;k++){
    if(abs(data[k]-(coordX[k]*coordY[k]*coordY[k]+1.0))>0.01){
      PRINT("ERROR" << "(" << rankX << "," << rankY <<","<<rank<< "), k:" << k<<" received: "<<data[k]<<" expected:"<<coordX[k]*coordY[k]*coordY[k]+1.0);
      hasFailed=true;
    }
  }
  if(hasFailed){
	PRINT("FAIL!");
  }else{ 
  	PRINT("SUCCESS!");
  }
  //destroy_(daemon_args,false);  
  
  MPI_Finalize();
  std::cout << "Stop communication proxy" << std::endl;
  return 0;
}

