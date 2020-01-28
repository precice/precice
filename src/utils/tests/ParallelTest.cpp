#include <string>
#include "testing/Testing.hpp"
#include "utils/Parallel.hpp"

using namespace precice;

BOOST_AUTO_TEST_SUITE(UtilsTests)

BOOST_AUTO_TEST_SUITE(Parallel)

#ifndef PRECICE_NO_MPI

BOOST_AUTO_TEST_CASE(RestrictCommTest)
{
  PRECICE_TEST(4_ranks);
  using Par     = utils::Parallel;

  auto baseComm = Par::current();
  Par::restrictCommunicator(3);
  auto restComm = Par::current();
  BOOST_TEST(restComm);
  BOOST_TEST(restComm->parent == baseComm);

  if (context.rank == 3) {
    BOOST_TEST(restComm->isNull());
  } else {
    BOOST_TEST(!restComm->isNull());
    BOOST_TEST(restComm->size() == 3);
    BOOST_TEST(restComm->rank() <= 2);
  }
}

BOOST_AUTO_TEST_CASE(SplitCommTest)
{
  PRECICE_TEST(3_ranks);
  using Par     = utils::Parallel;

  std::string name{"Group"};
  name.append(context.rank == 2 ? "Two" : "One");
  auto baseComm = Par::current();
  Par::splitCommunicator(name);
  auto splitComm = Par::current();
  BOOST_TEST(splitComm->parent == baseComm);

  BOOST_TEST(!splitComm->isNull());

  const auto &groups = splitComm->groups;
  BOOST_TEST(groups.size() == 2);
  BOOST_TEST(groups[0].id == 0);
  BOOST_TEST(groups[1].id == 1);
  BOOST_TEST(groups[0].name == std::string("GroupOne"));
  BOOST_TEST(groups[1].name == std::string("GroupTwo"));
  BOOST_TEST(groups[0].leaderRank == 0);
  BOOST_TEST(groups[1].leaderRank == 2);
  BOOST_TEST(groups[0].size == 2);
  BOOST_TEST(groups[1].size == 1);
}

#endif // not PRECICE_NO_MPI

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
