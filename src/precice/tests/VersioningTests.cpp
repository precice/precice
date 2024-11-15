#include <algorithm>
#include <boost/test/tools/old/interface.hpp>
#include <ostream>
#include <string>
#include <type_traits>

#include "precice/Tooling.hpp"
#include "precice/Version.h"
#include "testing/Testing.hpp"
#include "utils/assertion.hpp"

BOOST_AUTO_TEST_SUITE(PreciceTests)
BOOST_AUTO_TEST_SUITE(Versioning)

BOOST_AUTO_TEST_CASE(VersionInformation)
{
  PRECICE_TEST(1_rank);
  const std::string info{::precice::getVersionInformation()};
  BOOST_TEST_CONTEXT("Info is \"" << info << "\"")
  {
    BOOST_TEST(!info.empty(), "The info contains at least the version.");

    auto semis = std::count(info.begin(), info.end(), ';');
    BOOST_TEST(semis >= 2, "The info contains " << semis << " of at least 3 sections ");

    auto firstSemi = info.find(';');
    auto version   = info.substr(0, firstSemi);
    BOOST_TEST_INFO("Section: " << version);
    BOOST_TEST(!version.empty(), "The info contains the version section");

    auto secondSemi = info.find(';', firstSemi + 1);
    auto revision   = info.substr(firstSemi + 1, secondSemi - firstSemi - 1);
    BOOST_TEST_INFO("Section: " << revision);
    BOOST_TEST(!revision.empty(), "The info contains the revision section");
  }
}

BOOST_AUTO_TEST_CASE(VersionMacros)
{
  PRECICE_TEST(1_rank);

  BOOST_REQUIRE(std::is_integral<decltype(PRECICE_VERSION_MAJOR)>::value);
  BOOST_REQUIRE(std::is_integral<decltype(PRECICE_VERSION_MINOR)>::value);
  BOOST_REQUIRE(std::is_integral<decltype(PRECICE_VERSION_PATCH)>::value);

  BOOST_TEST(PRECICE_VERSION_MAJOR > 0);
  BOOST_TEST(PRECICE_VERSION_MINOR >= 0);
  BOOST_TEST(PRECICE_VERSION_PATCH >= 0);

  BOOST_REQUIRE_NO_THROW(std::string{PRECICE_VERSION});

  const std::string version{PRECICE_VERSION};
  // Minimal length of the version string is 5 chars X.Y.Z
  BOOST_REQUIRE(version.length() >= 5);
  // The version string contains 2 dots
  BOOST_REQUIRE(std::count(version.begin(), version.end(), '.') == 2);
}

BOOST_AUTO_TEST_CASE(VersionMacroGreaterEqual)
{
  PRECICE_TEST(1_rank);

  BOOST_REQUIRE(PRECICE_VERSION_GREATER_EQUAL(0, 0, 0));
  BOOST_REQUIRE(PRECICE_VERSION_GREATER_EQUAL(1, 0, 0));
  BOOST_REQUIRE(PRECICE_VERSION_GREATER_EQUAL(1, 1, 1));

  // Equality
  BOOST_REQUIRE(PRECICE_VERSION_GREATER_EQUAL(
      PRECICE_VERSION_MAJOR,
      PRECICE_VERSION_MINOR,
      PRECICE_VERSION_PATCH));

  // Wild guess
  BOOST_REQUIRE(!PRECICE_VERSION_GREATER_EQUAL(99, 9, 99));
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
