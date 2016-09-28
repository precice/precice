#ifndef PRECICE_IO_TESTS_EXPORTCONFIGURATIONTEST_HPP_
#define PRECICE_IO_TESTS_EXPORTCONFIGURATIONTEST_HPP_

#include "tarch/tests/TestCase.h"
#include "logging/Logger.hpp"

namespace precice {
namespace io {
namespace tests {

class ExportConfigurationTest : public tarch::tests::TestCase
{
public:

   ExportConfigurationTest ();

   virtual ~ExportConfigurationTest() {};

   /**
    * @brief Empty.
    */
   virtual void setUp () {}

   virtual void run ();

private:

   static logging::Logger _log;

   void testConfiguration ();
};

}}} // namespace precice, io, tests

#endif /* PRECICE_IO_TESTS_EXPORTCONFIGURATIONTEST_HPP_ */
