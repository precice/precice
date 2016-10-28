#pragma once

#ifndef PRECICE_NO_MPI

#include "../impl/ParallelMatrixOperations.hpp"

#include "com/Communication.hpp"
#include "m2n/M2N.hpp"
#include "tarch/tests/TestCase.h"
#include "logging/Logger.hpp"
#include <vector>

namespace precice {
namespace cplscheme {
namespace tests {

/**
 * @brief Provides tests for class ParallelMatrixOperationsTest for master-slave mode.
 */
class ParallelMatrixOperationsTest : public tarch::tests::TestCase
{
public:

    /**
     * @brief Constructor.
     */
	ParallelMatrixOperationsTest ();

   /**
    * @brief Destructor.
    */
   virtual ~ParallelMatrixOperationsTest() {};

   /**
    * @brief Empty.
    */
   virtual void setUp () {}

   /**
    * @brief Calls all test methods.
    */
   virtual void run ();


private:

   static logging::Logger _log;

   /**
    * @brief Tests the correct postprocessing for MVQN-like vector data
    */
   void testParallelMatrixMatrixOp();

   void testParVectorOperations();

   void validate_result_equals_reference(
   			Eigen::MatrixXd& result_local,
   			Eigen::MatrixXd& reference_global,
   			std::vector<int>& offsets,
   			bool partitionedRowWise);   
};

}}} // namespace precice, cplscheme, tests

#endif // PRECICE_NO_MPI
