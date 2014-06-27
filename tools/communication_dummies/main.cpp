#include "mpi.h"
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>

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

int main(int argc, char** argv)
{
  std::cout << "Running communication proxy" << std::endl;
  MPI_Init(&argc, &argv);

  if (argc == 1 || argc != 4){
    PRINT("Usage: ./proxy proX proY first/second=0/1");
    return 1;
  }

  int proX = atoi(argv[1]);
  int proY = atoi(argv[2]);
  bool doesFirstStep = not atof(argv[3]);

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
      coordX[(j-1)*chunkX +i-1] = (i-1)*dx + offsetX;
      coordY[(j-1)*chunkX +i-1] = (j-1)*dy + offsetY;
    }
  }


  if(doesFirstStep){
    PRINT("First Step");

    for(int k=0;k<pointSize;k++){
      data[k]=coordX[k]*coordY[k]*coordY[k];
    }

    //send data

  }
  else{
    PRINT("Second Step");
    for(int k=0;k<pointSize;k++){
      data[k]=0.0;
    }

    //receive data

    for(int k=0;k<pointSize;k++){
      if(data[k]!=coordX[k]*coordY[k]*coordY[k]){
        PRINT("ERROR" << "(" << rankX << "," << rankY << "), k:" << k);
      }
    }

  }

//  std::ostringstream myStream;
//  myStream << "(" << rankX << "," << rankY << ")";
//  for(int k=0;k<pointSize;k++){
//    myStream << data[k] << ", ";
//  }
//  std::cout << myStream.str() << std::endl;




  //MPI_Wtime();


  MPI_Finalize();
  std::cout << "Stop communication proxy" << std::endl;
  return 0;
}

