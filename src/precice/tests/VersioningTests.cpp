#include <algorithm>
#include <ostream>
#include <string>
#include "precice/SolverInterface.hpp"
#include "testing/Testing.hpp"
#include "utils/assertion.hpp"

BOOST_AUTO_TEST_SUITE(PreciceTests)
BOOST_AUTO_TEST_SUITE(Versioning)

BOOST_AUTO_TEST_CASE(VersionInformation)
{
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

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
