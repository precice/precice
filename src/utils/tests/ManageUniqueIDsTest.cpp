#include "testing/Testing.hpp"
#include "utils/ManageUniqueIDs.hpp"

using namespace precice::utils;

BOOST_AUTO_TEST_SUITE(UtilsTests)

BOOST_AUTO_TEST_CASE(UniqueIDs)
{
  ManageUniqueIDs uniqueIDs;
  int             id = uniqueIDs.getFreeID();
  BOOST_TEST(id == 0);
  id = uniqueIDs.getFreeID();
  BOOST_TEST(id == 1);
  bool success = uniqueIDs.insertID(2);
  BOOST_TEST(success);
  id = uniqueIDs.getFreeID();
  BOOST_TEST(id == 3);
}

BOOST_AUTO_TEST_SUITE_END()
