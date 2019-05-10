#include "testing/Testing.hpp"
#include "testing/Fixtures.hpp"

#include "precice/SolverInterface.hpp"
#include "versions.hpp"
#include <sstream>

using namespace precice;


BOOST_AUTO_TEST_SUITE(PreciceTests)
BOOST_AUTO_TEST_SUITE(Miscellaneous)

BOOST_AUTO_TEST_CASE(TestGetVersion)
{
    auto version = SolverInterface::getVersion();
    BOOST_TEST(version.major >= 0);
    BOOST_TEST(version.minor >= 0);
    BOOST_TEST(version.patch >= 0);

    std::ostringstream oss;
    oss << version.major << '.'
        << version.minor << '.'
        << version.patch;
    BOOST_TEST(PRECICE_VERSION == oss.str());
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
