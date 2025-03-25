#include <array>
#include <string>
#include <string_view>
#include <vector>

#include "precice/span.hpp"
#include "testing/Testing.hpp"

BOOST_AUTO_TEST_SUITE(PreciceTests)
BOOST_AUTO_TEST_SUITE(Span)

using precice::span;

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(FromDoubleVector)
{
  PRECICE_TEST();
  std::vector<double> v{1.0, 2.0, 3.0, 4.0};
  span<const double>  const_span{v};
  BOOST_TEST(v == const_span, boost::test_tools::per_element());
  span<double> mut_span{v};
  BOOST_TEST(v == mut_span, boost::test_tools::per_element());
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(FromIntVector)
{
  PRECICE_TEST();
  std::vector<int> v{1, 2, 3, 4};
  span<const int>  const_span{v};
  BOOST_TEST(v == const_span, boost::test_tools::per_element());
  span<int> mut_span{v};
  BOOST_TEST(v == mut_span, boost::test_tools::per_element());
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(FromString)
{
  PRECICE_TEST();
  std::string      s = "hello there";
  span<const char> const_span{s};
  BOOST_TEST(s == std::string(const_span.data(), const_span.size()));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(FromCharVec)
{
  PRECICE_TEST();
  std::vector<char> s{'h', 'e', 'l', 'l', 'o', ' ', 't', 'h', 'e', 'r', 'e'};
  span<const char>  const_span{s};
  BOOST_CHECK_EQUAL_COLLECTIONS(
      s.begin(), s.end(),
      const_span.begin(), const_span.end());
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(FromLiteral)
{
  PRECICE_TEST();
  std::string      s = "hello there";
  span<const char> const_span{"hello there"}; // Includes the NULL byte

  BOOST_CHECK_EQUAL_COLLECTIONS(s.begin(), ++s.end(),
                                const_span.begin(), const_span.end());
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(FromStdArray)
{
  PRECICE_TEST();
  std::array<char, 11> s{'h', 'e', 'l', 'l', 'o', ' ', 't', 'h', 'e', 'r', 'e'};
  span<const char>     const_span{s};

  BOOST_CHECK_EQUAL_COLLECTIONS(
      s.begin(), s.end(),
      const_span.begin(), const_span.end());
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(FromCArray)
{
  PRECICE_TEST();
  char             s[] = "hello there";
  span<const char> const_span{s}; // Includes the NULL byte

  BOOST_CHECK_EQUAL_COLLECTIONS(
      s, s + 12,
      const_span.begin(), const_span.end());
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(FromStringRef)
{
  PRECICE_TEST();
  std::string      s  = "hello there";
  std::string     &sr = s;
  span<const char> const_span{sr};
  BOOST_TEST(s == std::string(const_span.data(), const_span.size()));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(FromConstStringRef)
{
  PRECICE_TEST();
  std::string        s   = "hello there";
  const std::string &scr = s;
  span<const char>   const_span{scr};
  BOOST_TEST(s == std::string(const_span.data(), const_span.size()));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(FromCString)
{
  PRECICE_TEST();
  std::string      s  = "hello there";
  char            *sp = s.data();
  span<const char> const_span{sp};
  BOOST_TEST(s == std::string(const_span.data(), const_span.size()));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(FromConstCString)
{
  PRECICE_TEST();
  const char      *s = "hello there";
  span<const char> const_span{s};
  BOOST_TEST(s == std::string(const_span.data(), const_span.size()));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(ToStringView)
{
  PRECICE_TEST();
  const char      *s = "hello there";
  span<const char> const_span{s};
  std::string_view sv{const_span.data(), const_span.size()};

  BOOST_TEST(s == std::string(const_span.data(), const_span.size()));
  BOOST_TEST(s == sv);
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
