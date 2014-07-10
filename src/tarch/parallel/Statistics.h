/*
 * h
 *
 *  Created on: 16.04.2013
 *      Author: liebm
 */

#ifndef STATISTICS_H_
#define STATISTICS_H_

#include "tarch/parallel/Node.h"
#include <mpi.h>
#include <climits>
#include <cfloat>
#include <cstdlib>

namespace tarch {
  namespace parallel {
    class Statistics;
  }}

class tarch::parallel::Statistics {
public:
  Statistics();

  void getAverageMinMax(
      int myValue,
      std::string& info
  ) {
    std::ostringstream msg;
    msg <<info;

  }

  /**
   * this method computes the average, maximum and minimum value
   * among all ranks.
   * this is the integer version.
   */
  static void addStatisticsToStream(
      int myValue,
      std::string unitInfo,
      std::ostringstream& msg
  ) {
    const int numberOfNodes = tarch::parallel::Node::getInstance().getNumberOfNodes();
    // receive buffer
    int* rbuf = 0;

    if(tarch::parallel::Node::getInstance().isGlobalMaster())
      rbuf = (int *)malloc(numberOfNodes*sizeof(int));

    MPI_Gather(&myValue, 1, MPI_INT, rbuf, 1, MPI_INT, 0, tarch::parallel::Node::getInstance().getCommunicator());

    if(tarch::parallel::Node::getInstance().isGlobalMaster()) {
      int max=0; int min=INT_MAX;
      int nodeMax = INT_MIN; int nodeMin = INT_MIN;
      long average = 0;
      for(int i=0; i<numberOfNodes; i++) {
        if(rbuf[i] > max) {
          max = rbuf[i];
          nodeMax = i;
        }
        if (rbuf[i] < min) {
          min = rbuf[i];
          nodeMin = i;
        }
        average += rbuf[i];
      }
      free(rbuf);
      average = average / numberOfNodes;
      msg <<" average: " <<average <<unitInfo <<", "
          "maximum: " <<max <<unitInfo <<" (Node " <<nodeMax <<") , "
          "minimum: " <<min <<unitInfo <<" (Node " <<nodeMin <<")";
    }
  }

  /**
   * this method computes the average, maximum and minimum value
   * among all ranks.
   * this is the double version.
   */
  static void addStatisticsToStream(
      double myValue,
      std::string unitInfo,
      std::ostringstream& msg
  ) {
    const int numberOfNodes = tarch::parallel::Node::getInstance().getNumberOfNodes();
    // receive buffer
    double* rbuf = 0;

    if(tarch::parallel::Node::getInstance().isGlobalMaster())
      rbuf = (double *)malloc(numberOfNodes*sizeof(double));

    MPI_Gather(&myValue, 1, MPI_DOUBLE, rbuf, 1, MPI_DOUBLE, 0, tarch::parallel::Node::getInstance().getCommunicator());

    if(tarch::parallel::Node::getInstance().isGlobalMaster()) {
      double min = DBL_MAX, max = -DBL_MAX, average=0;
      int nodeMin= INT_MIN, nodeMax= INT_MIN;

      for(int i=0; i<numberOfNodes; i++) {
        if(rbuf[i] > max) {
          max = rbuf[i];
          nodeMax = i;
        }
        if (rbuf[i] < min) {
          min = rbuf[i];
          nodeMin = i;
        }
        average += rbuf[i];
      }
      free(rbuf);
      average = average / numberOfNodes;
      msg <<" average: " <<average <<unitInfo <<", "
          "maximum: " <<max <<unitInfo <<" (Node " <<nodeMax <<") , "
          "minimum: " <<min <<unitInfo <<" (Node " <<nodeMin <<")";
    }
  }

  virtual
  ~Statistics();
};
#endif /* STATISTICS_H_ */
