#ifndef PRECICE_UTILS_PARALLELTEST_HPP_
#define PRECICE_UTILS_PARALLELTEST_HPP_

#include "tarch/tests/TestCase.h"
#include "tarch/logging/Log.h"

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

  static tarch::logging::Log _log;
};

}}} // namespace precice, utils, tests

#endif /* PARALLELTEST_HPP_ */
