#include <string>
#include <string_view>

#include "testing/Testing.hpp"
#include "utils/span.hpp"

BOOST_AUTO_TEST_SUITE(PreciceTests)
BOOST_AUTO_TEST_SUITE(Span)

using precice::span;

BOOST_AUTO_TEST_CASE(FromDoubleVector)
{
  PRECICE_TEST(1_rank);
  std::vector<double> v{1.0, 2.0, 3.0, 4.0};
  span<const double>  const_span{v};
  BOOST_TEST(v == const_span, boost::test_tools::per_element());
  span<double> mut_span{v};
  BOOST_TEST(v == mut_span, boost::test_tools::per_element());
}

BOOST_AUTO_TEST_CASE(FromIntVector)
{
  PRECICE_TEST(1_rank);
  std::vector<int> v{1, 2, 3, 4};
  span<const int>  const_span{v};
  BOOST_TEST(v == const_span, boost::test_tools::per_element());
  span<int> mut_span{v};
  BOOST_TEST(v == mut_span, boost::test_tools::per_element());
}

BOOST_AUTO_TEST_CASE(FromString)
{
  PRECICE_TEST(1_rank);
  std::string      s = "hello there";
  span<const char> const_span{s};
  BOOST_TEST(s == std::string(const_span.data(), const_span.size()));
}

BOOST_AUTO_TEST_CASE(FromCString)
{
  PRECICE_TEST(1_rank);
  const char *     s = "hello there";
  span<const char> const_span{s};
  BOOST_TEST(s == std::string(const_span.data(), const_span.size()));
}

BOOST_AUTO_TEST_CASE(ToStringView)
{
  PRECICE_TEST(1_rank);
  const char *     s = "hello there";
  span<const char> const_span{s};
  std::string_view sv{const_span.data(), const_span.size()};

  BOOST_TEST(s == std::string(const_span.data(), const_span.size()));
  BOOST_TEST(s == sv);
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
