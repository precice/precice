#pragma once

#include "serial/explicit/helpers.hpp"
#include "testing/Testing.hpp"

#include <boost/test/tools/detail/per_element_manip.hpp>
#include <precice/SolverInterface.hpp>
#include <vector>

namespace tests {
inline void cyclicExplicit(TestContext &context)
{
  precice::SolverInterface interface{
      context.name, context.config(), 0, 1};

  auto mid = "M" + context.name;
  auto rid = std::map<std::string, std::string>{{"A", "DCA"}, {"B", "DAB"}, {"C", "DBC"}}.at(context.name); //  mid
  auto wid = std::map<std::string, std::string>{{"A", "DAB"}, {"B", "DBC"}, {"C", "DCA"}}.at(context.name); //  mid

  // create mesh
  const std::vector<double> coords{1, 0, 2, 0};

  std::vector<int> ids(2);
  interface.setMeshVertices(mid, 2, coords.data(), ids.data());

  std::vector<double> data{0, 0};

  interface.writeBlockScalarData(wid, 2, ids.data(), data.data());
  interface.initialize();

  BOOST_TEST_CONTEXT("first step")
  {
    BOOST_TEST_REQUIRE(interface.isCouplingOngoing());
    interface.readBlockScalarData(rid, 2, ids.data(), data.data());
    std::vector<double> expected(2, (context.isNamed("A") ? 0 : 1));
    BOOST_TEST(data == expected, boost::test_tools::per_element());

    std::fill(data.begin(), data.end(), 1);
    expected = data;
    interface.writeBlockScalarData(wid, 2, ids.data(), data.data());
    interface.advance(1);
  }

  BOOST_TEST_CONTEXT("second step")
  {
    BOOST_TEST_REQUIRE(interface.isCouplingOngoing());
    interface.readBlockScalarData(rid, 2, ids.data(), data.data());
    std::vector<double> expected(2, (context.isNamed("A") ? 1 : 2));
    BOOST_TEST(data == expected, boost::test_tools::per_element());

    std::fill(data.begin(), data.end(), 2);
    interface.writeBlockScalarData(wid, 2, ids.data(), data.data());
    interface.advance(1);
  }

  BOOST_TEST_REQUIRE(!interface.isCouplingOngoing());
}
} // namespace tests
