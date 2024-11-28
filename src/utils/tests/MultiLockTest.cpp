#include <map>
#include <string>
#include "math/constants.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"
#include "utils/MultiLock.hpp"
#include "utils/String.hpp"

using namespace precice;
using namespace precice::utils;

BOOST_AUTO_TEST_SUITE(UtilsTests)
BOOST_AUTO_TEST_SUITE(MultiLockTests)

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(MultiLockTest)
{
  PRECICE_TEST();
  MultiLock<std::string> mlock;

  BOOST_TEST(mlock.checkAll());
  BOOST_TEST(!mlock.contains("A"));
  BOOST_CHECK_THROW(mlock.check("A"), LockNotFoundException);
  BOOST_CHECK_THROW(mlock.lock("A"), LockNotFoundException);
  BOOST_CHECK_THROW(mlock.unlock("A"), LockNotFoundException);
  BOOST_TEST(mlock.size() == 0);
  BOOST_TEST(mlock.countUnlocked() == 0);
  BOOST_TEST(mlock.countLocked() == 0);

  mlock.add("A", true);
  BOOST_TEST(mlock.contains("A"));
  BOOST_TEST(!mlock.contains("B"));
  BOOST_TEST(mlock.check("A"));
  BOOST_TEST(mlock.checkAll());
  BOOST_TEST(mlock.size() == 1);
  BOOST_TEST(mlock.countUnlocked() == 0);
  BOOST_TEST(mlock.countLocked() == 1);

  mlock.add("B", false);
  BOOST_TEST(mlock.check("A"));
  BOOST_TEST(mlock.contains("B"));
  BOOST_TEST(!mlock.check("B"));
  BOOST_TEST(!mlock.checkAll());
  BOOST_TEST(mlock.size() == 2);
  BOOST_TEST(mlock.countUnlocked() == 1);
  BOOST_TEST(mlock.countLocked() == 1);

  mlock.lock("B");
  BOOST_TEST(mlock.check("A"));
  BOOST_TEST(mlock.check("B"));
  BOOST_TEST(mlock.checkAll());
  BOOST_TEST(mlock.size() == 2);
  BOOST_TEST(mlock.countUnlocked() == 0);
  BOOST_TEST(mlock.countLocked() == 2);

  mlock.unlockAll();
  BOOST_TEST(!mlock.check("A"));
  BOOST_TEST(!mlock.check("B"));
  BOOST_TEST(!mlock.checkAll());
  BOOST_TEST(mlock.countUnlocked() == 2);
  BOOST_TEST(mlock.countLocked() == 0);

  mlock.lockAll();
  BOOST_TEST(mlock.check("A"));
  BOOST_TEST(mlock.check("B"));
  BOOST_TEST(mlock.checkAll());
  BOOST_TEST(mlock.countUnlocked() == 0);
  BOOST_TEST(mlock.countLocked() == 2);
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
