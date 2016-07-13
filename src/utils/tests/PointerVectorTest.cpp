#include "PointerVectorTest.hpp"
#include "../PointerVector.hpp"
#include "../Parallel.hpp"
#include "../Helpers.hpp"

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::utils::tests::PointerVectorTest)

namespace precice {
namespace utils {
namespace tests {

tarch::logging::Log PointerVectorTest:: _log ( "precice::utils::tests::PointerVectoTest" );

PointerVectorTest:: PointerVectorTest ()
:
   TestCase ( "utils::tests::PointerVectorTest" )
{}

void PointerVectorTest:: run ()
{
   PRECICE_MASTER_ONLY {
      preciceTrace ( "run()" );
      ptr_vector<double> ptrVector;

   }
}

}}} // namespace precice, utils, tests
