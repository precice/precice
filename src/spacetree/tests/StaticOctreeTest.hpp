// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_SPACETREE_STATICOCTREETEST_HPP_
#define PRECICE_SPACETREE_STATICOCTREETEST_HPP_

#include "tarch/tests/TestCase.h"
#include "tarch/logging/Log.h"
#include "spacetree/tests/SpacetreeTestScenarios.hpp"
#include "spacetree/StaticOctree.hpp"

namespace precice {
namespace spacetree {
namespace tests {

/**
 * @brief Provides test cases for class RegularSpacetree.
 */
class StaticOctreeTest : public tarch::tests::TestCase
{
public:

  /**
   * @brief Constructor.
   */
  //StaticOctreeTest();

  /**
   * @brief Empty.
   */
  virtual void setUp() {}

  /**
   * @brief Runs all tests.
   */
	virtual void run();

private:

	struct StaticOctreeFactory : public SpacetreeTestScenarios::SpacetreeFactory
	{
    virtual PtrSpacetree createSpacetree (
      const utils::DynVector& offset,
      const utils::DynVector& halflengths,
      double                  refinementLimit )
    {
      return PtrSpacetree(new StaticOctree(offset, halflengths[0], refinementLimit));
    }
	};

	// @brief Logging device.
	static tarch::logging::Log _log;
};

}}} // namespace precice, spacetree, tests

#endif /* PRECICE_SPACETREE_STATICOCTREETEST_HPP_ */
