#ifndef PRECICE_GEOMETRY_TESTS_SPHERETEST_HPP_
#define PRECICE_GEOMETRY_TESTS_SPHERETEST_HPP_

#include "tarch/tests/TestCase.h"
#include "logging/Logger.hpp"

namespace precice {
namespace geometry {
namespace tests {

class SphereTest : public tarch::tests::TestCase
{
public:

   static logging::Logger _log;

   SphereTest ();

   virtual ~SphereTest () {};

   /**
    * @brief Empty.
    */
   virtual void setUp () {}

   virtual void run ();
};

}}} // namespace precice, geometry, tests

#endif // PRECICE_GEOMETRY_TESTS_SPHERETEST_HPP_
