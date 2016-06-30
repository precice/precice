#ifndef PRECICE_IO_TESTS_TXTWRITERREADERTEST_HPP_
#define PRECICE_IO_TESTS_TXTWRITERREADERTEST_HPP_

#include "tarch/tests/TestCase.h"
#include "tarch/logging/Log.h"

namespace precice {
namespace io {
namespace tests {

class TXTWriterReaderTest : public tarch::tests::TestCase
{
public:

  /**
   * @brief Constructor.
   */
  TXTWriterReaderTest();

  /**
   * @brief Destructor, empty.
   */
  virtual ~TXTWriterReaderTest() {}

  /**
   * @brief Empty.
   */
  virtual void setUp() {}

  /**
   * @brief Calls all tests.
   */
  virtual void run();

private:

  static tarch::logging::Log _log;

  /**
   * @brief Tests writing table data to a file.
   */
  void test();
};

}}} // namespace precice, io, tests

#endif // PRECICE_IO_TESTS_TXTWRITERREADERTEST_HPP_
