#include "testing/Testing.hpp"
#include "utils/Hash.hpp"

BOOST_AUTO_TEST_SUITE(UtilsTests)
BOOST_AUTO_TEST_SUITE(HashTests)

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(Empty)
{
  PRECICE_TEST();
  BOOST_TEST(precice::utils::preciceHash("") == "4d572a0ba37b5812b95a6e7b87898eea");
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(Letter)
{
  PRECICE_TEST();
  BOOST_TEST(precice::utils::preciceHash("A") == "28c1219289105e64885266336db22f68");
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(Digit)
{
  PRECICE_TEST();
  BOOST_TEST(precice::utils::preciceHash("1") == "cecd78530a875345828b93a514cb61c1");
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(XMLTags)
{
  PRECICE_TEST();
  BOOST_TEST(precice::utils::preciceHash("<precice-configuration></precice-configuration>") == "ea6c463e64ba52f7b887e18452156fe6");
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
