#include <algorithm>
#include <boost/test/framework.hpp>
#include <cstdlib>
#include <string>

#include "logging/LogMacros.hpp"
#include "logging/Logger.hpp"
#include "testing/Testing.hpp"
#include "utils/assertion.hpp"

namespace precice {
namespace testing {

std::string getPathToRepository()
{
  precice::logging::Logger _log("testing");
  char *                   preciceRoot = std::getenv("PRECICE_ROOT");
  PRECICE_CHECK(preciceRoot != nullptr,
                "Environment variable PRECICE_ROOT is required to run the tests, but has not been set. Please set it to the root directory of the precice repository.");
  return std::string(preciceRoot);
}

std::string getPathToSources()
{
  return getPathToRepository() + "/src";
}

std::string getPathToTests()
{
  return getPathToRepository() + "/tests";
}

std::vector<std::string> getTestSuites()
{
  std::vector<std::string> suites;

  auto &current = boost::unit_test::framework::current_test_case();

  auto pid = current.p_parent_id;
  while (true) {
    auto &parent = boost::unit_test::framework::get<boost::unit_test::test_suite>(pid);
    // the parent type may be "case", "suite" or "module"
    // "case" will not appear as it can only be a leaf (current)
    // "module" is the root of the tree
    PRECICE_ASSERT(parent.p_type_name != "case");
    if (parent.p_type_name == "suite") {
      suites.emplace_back(parent.p_name);
      pid = parent.p_parent_id;
    } else {
      // We hit the root of the tree/"module".
      // The suites are in backwards order: leaf-root
      std::reverse(suites.begin(), suites.end());
      return suites;
    }
  }
}

std::string getTestName()
{
  return boost::unit_test::framework::current_test_case().p_name;
}

} // namespace testing
} // namespace precice
