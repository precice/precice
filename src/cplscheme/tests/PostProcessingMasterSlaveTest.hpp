#ifndef PRECICE_NO_MPI

#ifndef PRECICE_POSTPROCESSING_TESTS_MASTERSLAVEPOSTPROCESSING_HPP_
#define PRECICE_POSTPROCESSING_TESTS_MASTERSLAVEPOSTPROCESSING_HPP_

#include "com/Communication.hpp"
#include "m2n/M2N.hpp"
#include "tarch/tests/TestCase.h"
#include "tarch/logging/Log.h"
#include "utils/xml/XMLTag.hpp"
#include "cplscheme/SharedPointer.hpp"
#include <string>
#include <vector>

namespace precice {
namespace cplscheme {
namespace tests {

/**
 * @brief Provides tests for class IQNILSPostProcessing for master-slave mode.
 */
class PostProcessingMasterSlaveTest : public tarch::tests::TestCase
{
public:

    /**
     * @brief Constructor.
     */
	PostProcessingMasterSlaveTest ();

   /**
    * @brief Destructor.
    */
   virtual ~PostProcessingMasterSlaveTest() {};

   /**
    * @brief Empty.
    */
   virtual void setUp () {}

   /**
    * @brief Calls all test methods.
    */
   virtual void run ();

   typedef std::map<int,PtrCouplingData> DataMap;
   typedef tarch::la::DynamicColumnMatrix<double> DataMatrix;

private:

   static tarch::logging::Log _log;

# ifndef PRECICE_NO_MPI

   /**
    * @brief Tests the correct postprocessing for VIQN-like vector data
    */
   void testVIQNILSpp();

   /**
    * @brief Tests the correct postprocessing for MVQN-like vector data
    */
   void testVIQNIMVJpp();

   void testIMVJ_effUpdate_pp();

# endif // not PRECICE_NO_MPI
};

}}} // namespace precice, cplscheme, tests

#endif /* PRECICE_POSTPROCESSING_TESTS_MASTERSLAVEPOSTPROCESSING_HPP_ */
#endif // PRECICE_NO_MPI
