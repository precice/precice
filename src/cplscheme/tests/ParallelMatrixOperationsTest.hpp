/*
 * ParallelMatrixOperationsTest.hpp
 *
 *  Created on: Aug 24, 2015
 *      Author: scheufks
 */
#ifndef PRECICE_NO_MPI
#ifndef PARALLELMATRIXOPERATIONSTEST_HPP_
#define PARALLELMATRIXOPERATIONSTEST_HPP_

#include "../impl/ParallelMatrixOperations.hpp"

#include "com/Communication.hpp"
#include "m2n/M2N.hpp"
#include "tarch/tests/TestCase.h"
#include "logging/Logger.hpp"
#include <string>
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
   void testParallelMatrixMatrixOp_tarch();
   void testParallelMatrixMatrixOp_Eigen();

   void testParVectorOperations();

   void validate_result_equals_reference(
			tarch::la::DynamicMatrix<double>& result_local,
			tarch::la::DynamicMatrix<double>& reference_global,
			std::vector<int>& offsets,
			bool partitionedRowWise);
   void validate_result_equals_reference(
   			Eigen::MatrixXd& result_local,
   			Eigen::MatrixXd& reference_global,
   			std::vector<int>& offsets,
   			bool partitionedRowWise);
   void validate_result_equals_reference(
			tarch::la::DynamicVector<double>& result_local,
			tarch::la::DynamicVector<double>& reference_global,
			std::vector<int>& offsets);

};

}}} // namespace precice, cplscheme, tests


#endif /* PARALLELMATRIXOPERATIONSTEST_HPP_ */
#endif // PRECICE_NO_MPI
