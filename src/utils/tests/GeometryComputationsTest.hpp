// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_UTILS_TESTS_GEOMETRYCOMPUTATIONSTEST_HPP_
#define PRECICE_UTILS_TESTS_GEOMETRYCOMPUTATIONSTEST_HPP_

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

	void testCollinear();

	void testBetween();

	void testSegmentsIntersect();

	void testSegmentPlaneIntersection();

	void testProjectVector();

	void testContainedInTriangle();

	void testContainedInHyperrectangle();
};

}}} // namespace precice, utils, tests

#endif /* PRECICE_UTILS_TESTS_GEOMETRYCOMPUTATIONSTEST_HPP_ */
