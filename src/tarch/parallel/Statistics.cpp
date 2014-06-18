#ifndef PRECICE_NO_MPI
#include "mpi.h"
#endif
#include "Statistics.h"

//using namespace tarch::parallel;
//
//
//void addStatisticsToStream(
//    int myValue,
//    std::string unitInfo,
//    std::ostringstream& msg
//) {
//  const int numberOfNodes = Node::getInstance().getNumberOfNodes();
//  // receive buffer
//  int* rbuf = 0;
//
//  if(Node::getInstance().isGlobalMaster())
//    rbuf = (int *)malloc(numberOfNodes*sizeof(int));
//
//  MPI_Gather(&myValue, 1, MPI_INT, rbuf, 1, MPI_INT, 0, Node::getInstance().getCommunicator());
//
//  if(Node::getInstance().isGlobalMaster()) {
//    int max=0; int min=INT_MAX;
//    int nodeMax = INT_MIN; int nodeMin = INT_MIN;
//    long average = 0;
//    for(int i=0; i<numberOfNodes; i++) {
//      if(rbuf[i] > max) {
//        max = rbuf[i];
//        nodeMax = i;
//      }
//      if (rbuf[i] < min) {
//        min = rbuf[i];
//        nodeMin = i;
//      }
//      average += rbuf[i];
//    }
//    free(rbuf);
//    average = average / numberOfNodes;
//    msg <<" average: " <<average <<unitInfo <<", "
//        "maximum: " <<max <<unitInfo <<" (Node " <<nodeMax <<") , "
//        "minimum: " <<min <<unitInfo <<" (Node " <<nodeMin <<")";
//  }
//}
//
//void addStatisticsToStream(
//    double myValue,
//    std::string unitInfo,
//    std::ostringstream& msg
//) {
//  const int numberOfNodes = Node::getInstance().getNumberOfNodes();
//  // receive buffer
//  double* rbuf = 0;
//
//  if(Node::getInstance().isGlobalMaster())
//    rbuf = (double *)malloc(numberOfNodes*sizeof(double));
//
//  MPI_Gather(&myValue, 1, MPI_DOUBLE, rbuf, 1, MPI_DOUBLE, 0, Node::getInstance().getCommunicator());
//
//  if(Node::getInstance().isGlobalMaster()) {
//    double min = DBL_MAX, max = DBL_MIN, average=0;
//    int nodeMin= INT_MIN, nodeMax= INT_MIN;
//
//    for(int i=0; i<numberOfNodes; i++) {
//      if(rbuf[i] > max) {
//        max = rbuf[i];
//        nodeMax = i;
//      }
//      if (rbuf[i] < min) {
//        min = rbuf[i];
//        nodeMin = i;
//      }
//      average += rbuf[i];
//    }
//    free(rbuf);
//    average = average / numberOfNodes;
//    msg <<" average: " <<average <<unitInfo <<", "
//        "maximum: " <<max <<unitInfo <<" (Node " <<nodeMax <<") , "
//        "minimum: " <<min <<unitInfo <<" (Node " <<nodeMin <<")";
//  }
//}
