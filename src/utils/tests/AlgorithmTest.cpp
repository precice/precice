#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <functional>
#include <ostream>
#include <string>
#include <utility>
#include <vector>
#include "math/constants.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"
#include "utils/String.hpp"
#include "utils/algorithm.hpp"

using namespace precice;
namespace pu = precice::utils;

BOOST_AUTO_TEST_SUITE(UtilsTests)
BOOST_AUTO_TEST_SUITE(AlgorithmTests)

BOOST_AUTO_TEST_CASE(MakeArray)
{
  PRECICE_TEST(1_rank);
  auto a = pu::make_array(1, 2, 3);
  BOOST_TEST(a.size() == 3);
  BOOST_TEST(a.at(0) == 1);
  BOOST_TEST(a.at(1) == 2);
  BOOST_TEST(a.at(2) == 3);
}

BOOST_AUTO_TEST_CASE(UniqueElements)
{
  PRECICE_TEST(1_rank);
  std::vector<int> y{1, 2, 3, 4, 5, 6, 7, 8, 9};
  BOOST_TEST(pu::unique_elements(y));
  BOOST_TEST(pu::unique_elements(y, [](int l, int r) { return l == r; }));

  std::vector<int> n1{1, 2, 3, 4, 5, 6, 7, 8, 9, 9};
  BOOST_TEST(!pu::unique_elements(n1));
  BOOST_TEST(!pu::unique_elements(n1, [](int l, int r) { return l == r; }));

  std::vector<int> n2{1, 1, 3, 4, 5, 6, 7, 8, 9};
  BOOST_TEST(!pu::unique_elements(n2));
  BOOST_TEST(!pu::unique_elements(n2, [](int l, int r) { return l == r; }));

  std::vector<int> e;
  BOOST_TEST(pu::unique_elements(e));
  BOOST_TEST(pu::unique_elements(e, [](int l, int r) { return l == r; }));
}

BOOST_AUTO_TEST_CASE(UniqueEigenElements)
{
  PRECICE_TEST(1_rank);
  Eigen::VectorXd v1(3);
  v1 << 1.0, 0.1, 0.2;
  Eigen::VectorXd v2(3);
  v2 << 0.1, 1.0, 0.2;
  Eigen::VectorXd v3(3);
  v3 << 0.1, 0.2, 1.0;

  std::vector<Eigen::VectorXd> case1{};
  BOOST_TEST(pu::unique_elements(case1));

  std::vector<Eigen::VectorXd> case2{v1};
  BOOST_TEST(pu::unique_elements(case2));

  std::vector<Eigen::VectorXd> case3{v1, v2};
  BOOST_TEST(pu::unique_elements(case3));

  std::vector<Eigen::VectorXd> case4{v1, v2, v1};
  BOOST_TEST(!pu::unique_elements(case4));

  std::vector<Eigen::VectorXd> case5{v1, v2, v3};
  BOOST_TEST(pu::unique_elements(case5));

  std::vector<Eigen::VectorXd> case6{v1, v2, v3, v1};
  BOOST_TEST(!pu::unique_elements(case6));
}

BOOST_AUTO_TEST_CASE(Mismatch)
{
  PRECICE_TEST(1_rank);
  std::vector<int> a{1, 2, 3, 4, 5, 6, 7, 8, 9};
  std::vector<int> b{1, 2, 3, 4, 5, 0, 9};

  auto aa = pu::mismatch(
      a.begin(), a.end(),
      a.begin(), a.end());
  BOOST_TEST((aa.first == aa.second));
  BOOST_TEST((aa.first == a.end()));

  auto ab = pu::mismatch(
      a.begin(), a.end(),
      b.begin(), b.end());
  BOOST_TEST((ab.first != a.end()));
  BOOST_TEST(*ab.first == 6);
  BOOST_TEST(*ab.second == 0);
}

BOOST_AUTO_TEST_SUITE(RangePreview)

BOOST_AUTO_TEST_CASE(NormalRangePreview)
{
  PRECICE_TEST(1_rank);
  std::vector<int>   a{1, 2, 3, 4, 5, 6, 0};
  std::ostringstream oss;
  oss << pu::previewRange(2, a);
  std::string str{oss.str()};
  BOOST_TEST(str == "[1, 2, ... , 6, 0] min:0 max:6");
}

BOOST_AUTO_TEST_CASE(PrintNoElements)
{
  PRECICE_TEST(1_rank);
  std::vector<int>   a{1, 2, 3, 4, 5, 6, 0};
  std::ostringstream oss;
  oss << pu::previewRange(0, a);
  std::string str{oss.str()};
  BOOST_TEST(str == "[ ... ] min:0 max:6");
}

BOOST_AUTO_TEST_CASE(EmptyRange)
{
  PRECICE_TEST(1_rank);
  std::vector<int>   a;
  std::ostringstream oss;
  oss << pu::previewRange(3, a);
  std::string str{oss.str()};
  BOOST_TEST(str == "<Empty Range>");
}

BOOST_AUTO_TEST_SUITE_END() // Range

BOOST_AUTO_TEST_SUITE(ReorderArray)

BOOST_AUTO_TEST_CASE(OneElement)
{
  PRECICE_TEST(1_rank);
  std::array<int, 1> input{1};
  std::array<int, 1> order{0};
  std::array<int, 1> expected{1};

  auto reordered = utils::reorder_array(order, input);
  BOOST_TEST(reordered == expected);
}

BOOST_AUTO_TEST_CASE(AlreadySorted)
{
  PRECICE_TEST(1_rank);
  std::array<int, 3> input{3, 4, 5};
  std::array<int, 3> order{0, 1, 2};
  std::array<int, 3> expected{3, 4, 5};

  auto reordered = utils::reorder_array(order, input);
  BOOST_TEST(reordered == expected);
}

BOOST_AUTO_TEST_CASE(Reverse)
{
  PRECICE_TEST(1_rank);
  std::array<int, 3> input{3, 4, 5};
  std::array<int, 3> order{2, 1, 0};
  std::array<int, 3> expected{5, 4, 3};

  auto reordered = utils::reorder_array(order, input);
  BOOST_TEST(reordered == expected);
}

BOOST_AUTO_TEST_CASE(Scramble)
{
  PRECICE_TEST(1_rank);
  std::array<int, 3> input{3, 4, 5};
  std::array<int, 3> order{2, 0, 1};
  std::array<int, 3> expected{5, 3, 4};

  auto reordered = utils::reorder_array(order, input);
  BOOST_TEST(reordered == expected);
}

BOOST_AUTO_TEST_CASE(ScramblePointer)
{
  PRECICE_TEST(1_rank);
  int a = 1, b = 2;

  std::array<int *, 3> input{&a, &b, nullptr};
  std::array<int, 3>   order{2, 0, 1};
  std::array<int *, 3> expected{nullptr, &a, &b};

  auto reordered = utils::reorder_array(order, input);
  BOOST_TEST(reordered == expected);
}

BOOST_AUTO_TEST_SUITE_END() // ReorderArray

BOOST_AUTO_TEST_SUITE_END() // Alorithm

BOOST_AUTO_TEST_SUITE_END() // Utils
