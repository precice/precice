#ifndef PRECICE_IO_TESTS_TXTTABLEWRITER_HPP_
#define PRECICE_IO_TESTS_TXTTABLEWRITER_HPP_

#include "tarch/tests/TestCase.h"
#include "logging/Logger.hpp"

namespace precice {
namespace io {
namespace tests {

class TXTTableWriterTest : public tarch::tests::TestCase
{
public:

  /**
   * @brief Constructor.
   */
  TXTTableWriterTest ();

  /**
   * @brief Destructor, empty.
   */
  virtual ~TXTTableWriterTest () {}

  /**
   * @brief Empty.
   */
  virtual void setUp () {}

  /**
   * @brief Calls all tests.
   */
  virtual void run ();

private:

  static logging::Logger _log;

  /**
   * @brief Tests writing table data to a file.
   */
  void test ();
};

}}} // namespace precice, io, tests

#endif // PRECICE_IO_TESTS_TXTTABLEWRITER_HPP_
