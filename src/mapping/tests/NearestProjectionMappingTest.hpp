#ifndef PRECICE_MAPPING_NEARESTPROJETIONMAPPINGTEST_HPP_
#define PRECICE_MAPPING_NEARESTPROJETIONMAPPINGTEST_HPP_

#include "tarch/tests/TestCase.h"
#include "tarch/logging/Log.h"

namespace precice {
namespace mapping {
namespace tests {


class NearestProjectionMappingTest
:
   public tarch::tests::TestCase
{
public:

     NearestProjectionMappingTest ();

   virtual ~NearestProjectionMappingTest() {};

   virtual void setUp () {}

   virtual void run ();

private:

   static tarch::logging::Log _log;

   void testConservativeNonIncremental();

   void testConsistentNonIncremental();
};

}}} // namespace precice, mapping, tests

#endif /* PRECICE_MAPPING_NEARESTPROJETIONMAPPINGTEST_HPP_ */
