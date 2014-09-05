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

  //do we really need this here? might be a problem for solvers that create normal MPI?
  MPI_Init_thread(&argc,&argv,MPI_THREAD_MULTIPLE ,&provided);
  std::cout <<  "TM: " << MPI_THREAD_MULTIPLE << std::endl;


  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (size!=5){
    std::cout << "Please run with 5 mpi processes" << std::endl;
    return 1;
  }

  SocketAscodtCommunication com;

  if(rank==0){
    int vertexTable[] = {1,0,1,1,0,2,2,4,4,4};
    //maybe change to strings, compare comment in MainA.cpp
    int adressTable[] = {0,1,2,3,4};
    com.acceptConnection(vertexTable, adressTable);
  }

  int numberOfVertices;
  if(rank==0){
    numberOfVertices=2;
  }
  else if(rank==1){
    numberOfVertices=3;
  }
  else if(rank==2){
    numberOfVertices=2;
  }
  else if(rank==3){
    //rank 3 has no vertices
    numberOfVertices=0;
  }
  else if(rank==4){
    numberOfVertices=3;
  }

  //receive
  double *pData = NULL;

  if(rank==0){
    double data[] ={0.0,0.0};
    pData = &data[0];
    com.receive(data,2);
  }
  else if(rank==1){
    double data[] ={0.0,0.0,0.0};
    pData = &data[0];
    com.receive(data,3);
  }
  else if(rank==2){
    double data[] ={0.0,0.0};
    pData = &data[0];
    com.receive(data,2);
  }
  else if(rank==3){
    double data[] ={};
    pData = &data[0];
    com.receive(data,0);
  }
  else if(rank==4){
    double data[] ={0.0,0.0,0.0};
    pData = &data[0];
    com.receive(data,3);
  }
  //all receives are blocking, so all we have an implicit barrier here


  //TODO as an intermediate step, you could also add another test here,
  // if B already received the right data

  //modify data
  for(int i=0;i<numberOfVertices;i++){
    pData[i] = pData[i] + rank + 1;
  }

  //send data
  if(rank==0){
    com.send(pData,2);
  }
  else if(rank==1){
    com.send(pData,3);
  }
  else if(rank==2){
    com.send(pData,2);
  }
  else if(rank==3){
    com.send(pData,0);
  }
  else if(rank==4){
    com.send(pData,3);
  }

  MPI_Finalize();
  std::cout << "Stop communication dummy, rank: " << rank << std::endl;
  return 0;
}

