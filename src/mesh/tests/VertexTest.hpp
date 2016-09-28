#ifndef PRECICE_MESH_TESTS_VERTEXTEST_HPP_
#define PRECICE_MESH_TESTS_VERTEXTEST_HPP_

#include "tarch/tests/TestCase.h"
#include "logging/Logger.hpp"

namespace precice {
namespace mesh {
namespace tests {

/**
 * @brief Provides tests for class Vertex.
 */
class VertexTest : public tarch::tests::TestCase
{
public:

  /**
   * @brief Constructor.
   */
  VertexTest ();

  /**
   * @brief Destructor, empty.
   */
  virtual ~VertexTest () {};

  /**
   * @brief Setup for all tests, empty.
   */
  virtual void setUp () {}

  /**
   * @brief Runs all tests.
   */
  virtual void run ();

private:

  // @brief Logging device.
  static logging::Logger _log;

  /**
   * @brief Tests vertex.
   */
  void test ();
};

}}} // namespace precice, mesh, tests

#endif // PRECICE_MESH_TESTS_VERTEXTEST_HPP_
