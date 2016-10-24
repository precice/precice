#ifndef PRECICE_SPACETREE_DYNAMICOCTREETEST_HPP_
#define PRECICE_SPACETREE_DYNAMICOCTREETEST_HPP_

#include "tarch/tests/TestCase.h"
#include "logging/Logger.hpp"
#include "spacetree/DynamicOctree.hpp"
#include "spacetree/tests/SpacetreeTestScenarios.hpp"

namespace precice {
namespace spacetree {
namespace tests {

/**
 * @brief Provides test cases for class RegularSpacetree.
 */
class DynamicOctreeTest : public tarch::tests::TestCase
{
public:

  /**
   * @brief Constructor.
   */
  DynamicOctreeTest();

  /**
   * @brief Empty.
   */
  virtual void setUp() {}

  /**
   * @brief Runs all tests.
   */
  virtual void run();

private:

  struct DynamicOctreeFactory : public SpacetreeTestScenarios::SpacetreeFactory
  {
    virtual PtrSpacetree createSpacetree (
      const Eigen::VectorXd& offset,
      const Eigen::VectorXd& halflengths,
      double                 refinementLimit )
    {
      return PtrSpacetree(new DynamicOctree(offset, halflengths[0], refinementLimit));
    }
  };

  // @brief Logging device.
  static logging::Logger _log;


};

}}} // namespace precice, spacetree, tests

#endif /* PRECICE_SPACETREE_DYNAMICOCTREETEST_HPP_ */
