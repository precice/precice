#pragma once

#include "tarch/tests/TestCase.h"
#include "logging/Logger.hpp"

namespace precice {
namespace utils {
namespace tests {

/**
 * @brief Provides tests for methods in utils/GeometryComputations.
 */
class GeometryComputationsTest : public tarch::tests::TestCase
{
public:

  GeometryComputationsTest();

  virtual void setUp() {}

  virtual void run();

private:

  static logging::Logger _log;

  void testTriangleArea();
  
  void testTetraVolume();

  void testCollinear();

  void testBetween();

  void testSegmentsIntersect();

  void testSegmentPlaneIntersection();

  void testProjectVector();

  void testContainedInTriangle();

  void testContainedInHyperrectangle();
};

}}} // namespace precice, utils, tests

