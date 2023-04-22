#include <algorithm>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/test/framework.hpp>
#include <boost/test/tree/test_unit.hpp>
#include <cstdlib>
#include <string>

#include "logging/LogMacros.hpp"
#include "logging/Logger.hpp"
#include "testing/SourceLocation.hpp"
#include "testing/Testing.hpp"
#include "utils/assertion.hpp"

namespace precice::testing {

std::string getPathToRepository()
{
  precice::logging::Logger _log("testing");
  return precice::testing::SourceLocation;
}

std::string getPathToSources()
{
  return getPathToRepository() + "/src";
}

std::string getPathToTests()
{
  return getPathToRepository() + "/tests";
}

std::string getTestPath()
{
  const auto &cspan = boost::unit_test::framework::current_test_case().p_file_name;
  return {cspan.begin(), cspan.end()};
}

std::string getTestName()
{
  std::string name = boost::unit_test::framework::current_test_case().p_name;
  // check for data tests
  if (name[0] != '_') {
    return name;
  }
  auto parent = boost::unit_test::framework::current_test_case().p_parent_id;
  return boost::unit_test::framework::get<boost::unit_test::test_suite>(parent).p_name;
}

} // namespace precice::testing
