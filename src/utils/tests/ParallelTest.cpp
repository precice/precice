#include <string>
#include "testing/Testing.hpp"
#include "utils/Parallel.hpp"

using namespace precice;

BOOST_AUTO_TEST_SUITE(UtilsTests)

#ifndef PRECICE_NO_MPI

BOOST_AUTO_TEST_CASE(Parallel, *testing::MinRanks{3})
{
  using Par     = utils::Parallel;
  MPI_Comm comm = Par::getRestrictedCommunicator({0, 1, 2});
  if (Par::getProcessRank() <= 2) {
    Par::setGlobalCommunicator(comm);
    std::string group;
    int         rank = Par::getProcessRank();
    if ((rank == 0) || (rank == 1)) {
      group = "GroupOne";
    } else {
      BOOST_TEST(rank == 2, rank);
      group = "GroupTwo";
    }
    Par::splitCommunicator(group);

    const std::vector<Par::AccessorGroup> &groups = Par::getAccessorGroups();
    BOOST_TEST(groups.size() == 2);
    BOOST_TEST(groups[0].id == 0);
    BOOST_TEST(groups[1].id == 1);
    BOOST_TEST(groups[0].name == std::string("GroupOne"));
    BOOST_TEST(groups[1].name == std::string("GroupTwo"));
    BOOST_TEST(groups[0].leaderRank == 0);
    BOOST_TEST(groups[1].leaderRank == 2);
    BOOST_TEST(groups[0].size == 2);
    BOOST_TEST(groups[1].size == 1);

    Par::setGlobalCommunicator(Par::getCommunicatorWorld());
  }
}

#endif // not PRECICE_NO_MPI

BOOST_AUTO_TEST_SUITE_END()
