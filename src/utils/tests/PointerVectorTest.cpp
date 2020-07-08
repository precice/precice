#include "math/constants.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"
#include "utils/PointerVector.hpp"

using namespace precice;

BOOST_AUTO_TEST_SUITE(UtilsTests)

BOOST_AUTO_TEST_CASE(PointerVector)
{
  PRECICE_TEST(1_rank);
  utils::ptr_vector<double> ptrVector;
}

BOOST_AUTO_TEST_SUITE_END()
