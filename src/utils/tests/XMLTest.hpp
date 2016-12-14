#ifndef PRECICE_UTILS_TESTS_XMLTEST_HPP_
#define PRECICE_UTILS_TESTS_XMLTEST_HPP_

#include "tarch/tests/TestCase.h"
#include "logging/Logger.hpp"
#include "utils/xml/XMLTag.hpp"
#include "utils/Dimensions.hpp"

namespace precice {
namespace utils {
namespace tests {

/// Provides tests for classes in utils/xml/.
class XMLTest : public tarch::tests::TestCase, public XMLTag::Listener
{
public:

  XMLTest();

  /**
   * @brief Destructor, empty.
   */
  virtual ~XMLTest() {}

  /**
   * @brief Empty.
   */
  virtual void setUp();

  /**
   * @brief Runs all tests.
   */
  virtual void run();

  virtual void xmlTagCallback ( XMLTag& callingTag );

  virtual void xmlEndTagCallback ( utils::XMLTag& callingTag ) {}

private:

  static logging::Logger _log;

  std::string _testDirectory;

  utils::Vector2D _vector2D;
  utils::Vector3D _vector3D;
  utils::DynVector _dynVector;
  Eigen::VectorXd _eigenVectorXd;
  
  void testAttributeConcatenation();

  void testVectorAttributes();

  //void testNestedTags();
};

}}} // namespace precice, utils, tests

#endif /* PRECICE_UTILS_TESTS_XMLTEST_HPP_ */
