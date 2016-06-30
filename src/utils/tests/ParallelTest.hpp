#ifndef PRECICE_UTILS_PARALLELTEST_HPP_
#define PRECICE_UTILS_PARALLELTEST_HPP_

#include "tarch/tests/TestCase.h"
#include "logging/Logger.hpp"

namespace precice {
namespace utils {
namespace tests {

class ParallelTest : public tarch::tests::TestCase
{
public:

  ParallelTest ();

  virtual ~ParallelTest () {}

  virtual void setUp () {}

  virtual void run ();

private:

  static logging::Logger _log;
};

}}} // namespace precice, utils, tests

#endif /* PARALLELTEST_HPP_ */
