#include <boost/test/tools/context.hpp>
#include "testing/Testing.hpp"
#include "utils/MasterSlave.hpp"

using namespace precice;

BOOST_AUTO_TEST_SUITE(UtilsTests)

BOOST_AUTO_TEST_SUITE(MasterSlave)

BOOST_AUTO_TEST_CASE(SerialConfig)
{
  PRECICE_TEST(""_on(1_rank).setupMasterSlaves());

  BOOST_TEST(!utils::MasterSlave::isMaster());
  BOOST_TEST(!utils::MasterSlave::isSlave());
  BOOST_TEST(!utils::MasterSlave::isParallel());

  BOOST_TEST(utils::MasterSlave::getRank() == context.rank);
  BOOST_TEST(utils::MasterSlave::getSize() == context.size);

  { // ranks
    auto             ranksRange = utils::MasterSlave::allRanks();
    std::vector<int> ranks(ranksRange.begin(), ranksRange.end());
    BOOST_TEST(ranks.size() == 1);
    BOOST_TEST(ranks.front() == 0);
  }

  { // slaves
    auto slaves = utils::MasterSlave::allSlaves();
    BOOST_TEST((slaves.begin() == slaves.end()));
  }

  BOOST_TEST(!static_cast<bool>(utils::MasterSlave::_communication));
}

BOOST_AUTO_TEST_CASE(ParallelConfig)
{
  PRECICE_TEST(""_on(3_ranks).setupMasterSlaves());

  BOOST_TEST(utils::MasterSlave::isMaster() == context.isMaster());
  BOOST_TEST(utils::MasterSlave::isSlave() != context.isMaster());
  BOOST_TEST(utils::MasterSlave::isParallel());

  BOOST_TEST(utils::MasterSlave::getRank() == context.rank);
  BOOST_TEST(utils::MasterSlave::getSize() == context.size);

  { // ranks
    auto             ranksRange = utils::MasterSlave::allRanks();
    std::vector<int> ranks(ranksRange.begin(), ranksRange.end());
    std::vector<int> expected{0, 1, 2};
    BOOST_TEST(ranks == expected, boost::test_tools::per_element());
  }

  { // slaves
    auto             slaves = utils::MasterSlave::allSlaves();
    std::vector<int> ranks(slaves.begin(), slaves.end());
    std::vector<int> expected{1, 2};
    BOOST_TEST(ranks == expected, boost::test_tools::per_element());
  }

  BOOST_TEST(static_cast<bool>(utils::MasterSlave::_communication));
}

BOOST_AUTO_TEST_CASE(Parallell2norm)
{
  PRECICE_TEST(""_on(3_ranks).setupMasterSlaves());

  const double norm = 16.881943016134134;
  if (context.isMaster()) {
    Eigen::VectorXd v(3);
    v << 1, 2, 3;
    BOOST_TEST(utils::MasterSlave::l2norm(v) == norm);
  }
  if (context.isRank(1)) {
    Eigen::VectorXd v(2);
    v << 4, 5;
    BOOST_TEST(utils::MasterSlave::l2norm(v) == norm);
  }
  if (context.isRank(2)) {
    Eigen::VectorXd v(4);
    v << 6, 7, 8, 9;
    BOOST_TEST(utils::MasterSlave::l2norm(v) == norm);
  }
}

BOOST_AUTO_TEST_CASE(Seriall2norm)
{
  PRECICE_TEST(""_on(1_rank).setupMasterSlaves());

  const double    norm = 16.881943016134134;
  Eigen::VectorXd v(9);
  v << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  BOOST_TEST(utils::MasterSlave::l2norm(v) == norm);
}

BOOST_AUTO_TEST_CASE(Paralleldot)
{
  PRECICE_TEST(""_on(3_ranks).setupMasterSlaves());

  if (context.isMaster()) {
    Eigen::VectorXd u(3), v(3);
    u << 1, 2, 3;
    v << 9, 8, 7;
    BOOST_TEST(utils::MasterSlave::dot(u, v) == 165);
  }
  if (context.isRank(1)) {
    Eigen::VectorXd u(2), v(2);
    u << 4, 5;
    v << 6, 5;
    BOOST_TEST(utils::MasterSlave::dot(u, v) == 165);
  }
  if (context.isRank(2)) {
    Eigen::VectorXd u(4), v(4);
    u << 6, 7, 8, 9;
    v << 4, 3, 2, 1;
    BOOST_TEST(utils::MasterSlave::dot(u, v) == 165);
  }
}

BOOST_AUTO_TEST_CASE(Serialdot)
{
  PRECICE_TEST(""_on(1_rank).setupMasterSlaves());

  Eigen::VectorXd u(9), v(9);
  u << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  v << 9, 8, 7, 6, 5, 4, 3, 2, 1;
  BOOST_TEST(utils::MasterSlave::dot(u, v) == 165);
}

BOOST_AUTO_TEST_CASE(ParallelReduceSum)
{
  PRECICE_TEST(""_on(3_ranks).setupMasterSlaves());

  if (context.isMaster()) {
    {
      std::vector<double> in{1, 2, 3}, out{-1, -1, -1};
      utils::MasterSlave::reduceSum(in, out);
      BOOST_TEST(out == (std::vector<double>{12, 15, 18}), boost::test_tools::per_element());
    }
    {
      int in = 1, out = -1;
      utils::MasterSlave::reduceSum(in, out);
      BOOST_TEST(out == 9);
    }
    {
      double in = 1.1, out = -1;
      utils::MasterSlave::reduceSum(in, out);
      BOOST_TEST(out == 9.9);
    }
  }
  if (context.isRank(1)) {
    {
      std::vector<double> in{4, 5, 6}, out{-1, -1, -1};
      auto                expected = out;
      utils::MasterSlave::reduceSum(in, out);
      BOOST_TEST(out == expected, boost::test_tools::per_element());
    }
    {
      int in = 3, out = -1;
      utils::MasterSlave::reduceSum(in, out);
      BOOST_TEST(out == -1);
    }
    {
      double in = 3.3, out = -1;
      utils::MasterSlave::reduceSum(in, out);
      BOOST_TEST(out == -1);
    }
  }
  if (context.isRank(2)) {
    {
      std::vector<double> in{7, 8, 9}, out{-1, -1, -1};
      auto                expected = out;
      utils::MasterSlave::reduceSum(in, out);
      BOOST_TEST(out == expected, boost::test_tools::per_element());
    }
    {
      int in = 5, out = -1;
      utils::MasterSlave::reduceSum(in, out);
      BOOST_TEST(out == -1);
    }
    {
      double in = 5.5, out = -1;
      utils::MasterSlave::reduceSum(in, out);
      BOOST_TEST(out == -1);
    }
  }
}

BOOST_AUTO_TEST_CASE(SerialReduceSum)
{
  PRECICE_TEST(""_on(1_rank).setupMasterSlaves());

  {
    std::vector<double> in{1, 2, 3}, out{-1, -1, -1};
    utils::MasterSlave::reduceSum(in, out);
    BOOST_TEST(out == in, boost::test_tools::per_element());
  }
  {
    int in = 1, out = -1;
    utils::MasterSlave::reduceSum(in, out);
    BOOST_TEST(out == in);
  }
  {
    double in = 1.1, out = -1;
    utils::MasterSlave::reduceSum(in, out);
    BOOST_TEST(out == in);
  }
}

BOOST_AUTO_TEST_CASE(ParallelAllReduceSum)
{
  PRECICE_TEST(""_on(3_ranks).setupMasterSlaves());

  if (context.isMaster()) {
    {
      std::vector<double> in{1, 2, 3}, out{-1, -1, -1};
      utils::MasterSlave::allreduceSum(in, out);
      BOOST_TEST(out == (std::vector<double>{12, 15, 18}), boost::test_tools::per_element());
    }
    {
      int in = 1, out = -1;
      utils::MasterSlave::allreduceSum(in, out);
      BOOST_TEST(out == 9);
    }
    {
      double in = 1.1, out = -1;
      utils::MasterSlave::allreduceSum(in, out);
      BOOST_TEST(out == 9.9);
    }
  }
  if (context.isRank(1)) {
    {
      std::vector<double> in{4, 5, 6}, out{-1, -1, -1};
      utils::MasterSlave::allreduceSum(in, out);
      BOOST_TEST(out == (std::vector<double>{12, 15, 18}), boost::test_tools::per_element());
    }
    {
      int in = 3, out = -1;
      utils::MasterSlave::allreduceSum(in, out);
      BOOST_TEST(out == 9);
    }
    {
      double in = 3.3, out = -1;
      utils::MasterSlave::allreduceSum(in, out);
      BOOST_TEST(out == 9.9);
    }
  }
  if (context.isRank(2)) {
    {
      std::vector<double> in{7, 8, 9}, out{-1, -1, -1};
      utils::MasterSlave::allreduceSum(in, out);
      BOOST_TEST(out == (std::vector<double>{12, 15, 18}), boost::test_tools::per_element());
    }
    {
      int in = 5, out = -1;
      utils::MasterSlave::allreduceSum(in, out);
      BOOST_TEST(out == 9);
    }
    {
      double in = 5.5, out = -1;
      utils::MasterSlave::allreduceSum(in, out);
      BOOST_TEST(out == 9.9);
    }
  }
}

BOOST_AUTO_TEST_CASE(SerialAllReduceSum)
{
  PRECICE_TEST(""_on(1_rank).setupMasterSlaves());

  {
    std::vector<double> in{1, 2, 3}, out{-1, -1, -1};
    utils::MasterSlave::allreduceSum(in, out);
    BOOST_TEST(out == in, boost::test_tools::per_element());
  }
  {
    int in = 1, out = -1;
    utils::MasterSlave::allreduceSum(in, out);
    BOOST_TEST(out == in);
  }
  {
    double in = 1.1, out = -1;
    utils::MasterSlave::allreduceSum(in, out);
    BOOST_TEST(out == in);
  }
}

BOOST_AUTO_TEST_CASE(ParallelBroadcast)
{
  PRECICE_TEST(""_on(3_ranks).setupMasterSlaves());

  if (context.isMaster()) {
    {
      std::vector<double> in{1, 2, 3};
      utils::MasterSlave::broadcast(in);
      BOOST_TEST(in == (std::vector<double>{1, 2, 3}), boost::test_tools::per_element());
    }
    {
      bool in = true;
      utils::MasterSlave::broadcast(in);
      BOOST_TEST(in);
    }
    {
      double in = 9.9;
      utils::MasterSlave::broadcast(in);
      BOOST_TEST(in == 9.9);
    }
  } else {
    {
      std::vector<double> out{-1, -1, -1};
      utils::MasterSlave::broadcast(out);
      BOOST_TEST(out == (std::vector<double>{1, 2, 3}), boost::test_tools::per_element());
    }
    {
      bool out = false;
      utils::MasterSlave::broadcast(out);
      BOOST_TEST(out);
    }
    {
      double out = 9.9;
      utils::MasterSlave::broadcast(out);
      BOOST_TEST(out == 9.9);
    }
  }
}

BOOST_AUTO_TEST_CASE(SerialBroadcast)
{
  PRECICE_TEST(""_on(1_rank).setupMasterSlaves());
  {
    std::vector<double> in{1, 2, 3};
    utils::MasterSlave::broadcast(in);
    BOOST_TEST(in == (std::vector<double>{1, 2, 3}), boost::test_tools::per_element());
  }
  {
    bool in = true;
    utils::MasterSlave::broadcast(in);
    BOOST_TEST(in);
  }
  {
    double in = 9.9;
    utils::MasterSlave::broadcast(in);
    BOOST_TEST(in == 9.9);
  }
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
