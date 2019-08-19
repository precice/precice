#include <string>
#include "testing/Testing.hpp"
#include "utils/MultiLock.hpp"

using namespace precice::utils;

BOOST_AUTO_TEST_SUITE(UtilsTests)

BOOST_AUTO_TEST_CASE(MultiLockTest)
{
    MultiLock<std::string> mlock; 

    BOOST_TEST( mlock.checkAll());
    BOOST_TEST(!mlock.contains("A"));
    BOOST_CHECK_THROW(mlock.check("A"), LockNotFoundException);
    BOOST_CHECK_THROW(mlock.lock("A"), LockNotFoundException);
    BOOST_CHECK_THROW(mlock.unlock("A"), LockNotFoundException);

    mlock.add("A", true);
    BOOST_TEST( mlock.contains("A"));
    BOOST_TEST(!mlock.contains("B"));
    BOOST_TEST( mlock.check("A"));
    BOOST_TEST( mlock.checkAll());

    mlock.add("B", false);
    BOOST_TEST( mlock.check("A"));
    BOOST_TEST( mlock.contains("B"));
    BOOST_TEST(!mlock.check("B"));
    BOOST_TEST(!mlock.checkAll());

    mlock.lock("B");
    BOOST_TEST(mlock.check("A"));
    BOOST_TEST(mlock.check("B"));
    BOOST_TEST(mlock.checkAll());

    mlock.unlockAll();
    BOOST_TEST(!mlock.check("A"));
    BOOST_TEST(!mlock.check("B"));
    BOOST_TEST(!mlock.checkAll());

    mlock.lockAll();
    BOOST_TEST(mlock.check("A"));
    BOOST_TEST(mlock.check("B"));
    BOOST_TEST(mlock.checkAll());
}

BOOST_AUTO_TEST_SUITE_END()
