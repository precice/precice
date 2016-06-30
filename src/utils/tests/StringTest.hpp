#ifndef PRECICE_UTILS_STRINGTEST_HPP_
#define PRECICE_UTILS_STRINGTEST_HPP_

#include "tarch/tests/TestCase.h"
#include "logging/Logger.hpp"

namespace precice {
namespace utils {
namespace tests {

/**
 * @brief Provides tests for methods in utils/String.hpp.
 */
class StringTest : public tarch::tests::TestCase
{
public:

   StringTest();

   virtual ~StringTest() {}

   virtual void setUp() {}

   virtual void run();

private:

   static logging::Logger _log;

   void testWrapText();

   void testCheckAppendExtension();
};

}}} // namespace precice, utils, tests

#endif /* PRECICE_UTILS_STRINGTEST_HPP_ */
