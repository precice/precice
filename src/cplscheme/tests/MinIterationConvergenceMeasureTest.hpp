#ifndef PRECICE_CPLSCHEME_TEST_MINITERATIONCONVERGENCEMEASURETEST_HPP_
#define PRECICE_CPLSCHEME_TEST_MINITERATIONCONVERGENCEMEASURETEST_HPP_

#include "tarch/tests/TestCase.h"
#include "tarch/logging/Log.h"

namespace precice {
namespace cplscheme {
namespace tests {

class MinIterationConvergenceMeasureTest : public tarch::tests::TestCase
{
public:

   MinIterationConvergenceMeasureTest ();

   virtual ~MinIterationConvergenceMeasureTest () {};

   /**
    * @brief Empty.
    */
   virtual void setUp () {}

   virtual void run ();

private:

   static tarch::logging::Log _log;
};

}}} // namespace precice, cplscheme, tests

#endif // PRECICE_CPLSCHEME_TEST_MINITERATIONCONVERGENCEMEASURETEST_HPP_
