#include "mpi.h"
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>
#include "math.h"
#include "SocketAscodtCommunication.h"


int main(int argc, char** argv)
{
  std::cout << "Running communication dummy" << std::endl;
  int provided;

  //we should try in the end if this also works for the normal MPI,
  //then we would be on the safe side with the solver's MPI
  MPI_Init_thread(&argc,&argv,MPI_THREAD_MULTIPLE ,&provided);
  std::cout <<  "TM: " << MPI_THREAD_MULTIPLE << std::endl;


  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (size!=3){
    std::cout << "Please run with 3 mpi processes" << std::endl;
    return 1;
  }

  SocketAscodtCommunication com;

  if(rank==0){
    int vertexTable[] = {0,0,1,0,1,1,2,0,1,2};
    int adressTable[] = {0,1,2};
    // probably, we have to change adressTable to an array of strings, or 2 arrays,
    // whatever, feel free to adjust the interface here
    com.requestConnection(vertexTable, adressTable);
  }

  int numberOfVertices;
  if(rank==0){
    numberOfVertices=4;
  }
  else if(rank==1){
    numberOfVertices=4;
  }
  else if(rank==2){
    numberOfVertices=2;
  }

  // send and receive
  // both function should be blocking seen from here (not inside) -> use some sort of ack

  double *pData = NULL;

  if(rank==0){
    //data is sub-ordered just like vertexTable, so rank 0 has the 1st,
    //the 2nd, the 4th, and the 8th vertx
    double data[]={10.0,20.0,40.0,80.0};
    pData = &data[0];
    com.send(pData,4);
    com.receive(pData,4);
  }
  else if(rank==1){
    double data[]={30.0,50.0,60.0,90.0};
    pData = &data[0];
    com.send(pData,4);
    com.receive(pData,4);
  }
  else if(rank==2){
    double data[]={70.0,100.0};
    pData = &data[0];
    com.send(pData,2);
    com.receive(pData,2);
  }
  
  
  //check if everything worked out well
  bool correct = true;

  if(rank==0){
    double expectedData[]={12.0,21.0,42.0,85.0};
    for(int i=0;i<4;i++){
      if(pData[i]!=expectedData[i]){
        std::cout << "Error: " << pData[i] << " instead of " << expectedData[i]
                       << " on rank " << rank << std::endl;
        correct = false;
      }
    }
  }
  else if(rank==1){
    double expectedData[]={32.0,51.0,63.0,95.0};
    for(int i=0;i<4;i++){
      if(pData[i]!=expectedData[i]){
        std::cout << "Error: " << pData[i] << " instead of " << expectedData[i]
                       << " on rank " << rank << std::endl;
        correct = false;
      }
    }
  }
  else if(rank==2){
    double expectedData[]={73.0,105.0};
    for(int i=0;i<2;i++){
      if(pData[i]!=expectedData[i]){
        std::cout << "Error: " << pData[i] << " instead of " << expectedData[i]
                       << " on rank " << rank << std::endl;
        correct = false;
      }
    }
  }

  if(correct){
    std::cout << "Everything worked out well for rank: " << rank << std::endl;
  }

  MPI_Finalize();
  std::cout << "Stop communication dummy" << std::endl;
  return 0;
}

