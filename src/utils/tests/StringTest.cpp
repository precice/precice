#include <boost/test/tools/interface.hpp>
#include <boost/test/tools/old/interface.hpp>
#include <string>
#include "math/constants.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"
#include "utils/String.hpp"

using namespace precice;
using namespace precice::utils;

BOOST_AUTO_TEST_SUITE(UtilsTests)
BOOST_AUTO_TEST_SUITE(StringTests)

BOOST_AUTO_TEST_CASE(StringWrap)
{
  PRECICE_TEST(1_rank);
  std::string text("123456 1234567 12345678");
  std::string wrapped = wrapText(text, 4, 0);
  BOOST_TEST(wrapped == std::string("123456\n1234567\n12345678"));

  text    = "1234 1234 1234";
  wrapped = wrapText(text, 5, 0);
  BOOST_TEST(wrapped == std::string("1234\n1234\n1234"));

  text    = "1234 123 5";
  wrapped = wrapText(text, 5, 0);
  BOOST_TEST(wrapped == std::string("1234\n123 5"));

  text    = "1234 1234 1234";
  wrapped = wrapText(text, 7, 2);
  BOOST_TEST(wrapped == std::string("1234\n  1234\n  1234"));

  text    = "1234 1234 1234";
  wrapped = wrapText(text, 5, 2);
  BOOST_TEST(wrapped == std::string("1234\n  1234\n  1234"));

  text    = "12345678 1234 1 1234";
  wrapped = wrapText(text, 8, 2);
  BOOST_TEST(wrapped == std::string("12345678\n  1234 1\n  1234"));
}

BOOST_AUTO_TEST_CASE(StringAppendExtension)
{
  PRECICE_TEST(1_rank);
  std::string filename("somefile");
  std::string extension(".xyz");

  std::string result(filename);
  checkAppendExtension(result, extension);
  BOOST_TEST(result.compare(filename + extension) == 0);

  result = filename + extension;
  checkAppendExtension(result, extension);
  BOOST_TEST(result.compare(filename + extension) == 0);

  result = filename + extension + ".zyx";
  checkAppendExtension(result, extension);
  BOOST_TEST(result.compare(filename + extension + ".zyx" + extension) == 0);
}

BOOST_AUTO_TEST_CASE(ConvertStringToBool)
{
  PRECICE_TEST(1_rank);
  BOOST_TEST(convertStringToBool("tRUe") == true);
  BOOST_TEST(convertStringToBool("FALSE") == false);
  BOOST_TEST(convertStringToBool("oN") == true);
  BOOST_TEST(convertStringToBool("off") == false);
  BOOST_TEST(convertStringToBool("1") == true);
  BOOST_TEST(convertStringToBool("0") == false);
  BOOST_TEST(convertStringToBool("yes") == true);
  BOOST_TEST(convertStringToBool("no") == false);
}

BOOST_AUTO_TEST_CASE(StringMaker)
{
  PRECICE_TEST(1_rank);
  utils::StringMaker<10> sm;
  BOOST_REQUIRE(sm.data() != nullptr);
  BOOST_TEST(*sm.data() == '\0');
  BOOST_TEST(sm.str().empty());

  auto cstr = "1234";
  std::copy(cstr, cstr + 4, sm.data());
  BOOST_TEST(*sm.data() == '1');

  auto str = sm.str();
  BOOST_TEST(str.size() == 4);
  BOOST_TEST(str == cstr);

  sm.data()[2] = '\0';
  auto str2    = sm.str();
  BOOST_TEST(str2.size() == 2);
  BOOST_TEST(str2 == "12");
}

BOOST_AUTO_TEST_CASE(EditDistance)
{
  PRECICE_TEST(1_rank);
  using precice::utils::editDistance;
  BOOST_TEST(editDistance("left", "theft") == 2);
  BOOST_TEST(editDistance("aaa", "aa") == 1);
  BOOST_TEST(editDistance("aAa", "aaa") == 1);
  BOOST_TEST(editDistance("same", "same") == 0);
}

BOOST_AUTO_TEST_CASE(ComputeMatches)
{
  PRECICE_TEST(1_rank);
  std::set names{"left", "right", "up", "down"};
  auto     matches = utils::computeMatches("lewp", names);
  BOOST_TEST_REQUIRE(matches.size() == 4);

  BOOST_TEST(matches[0].name == "left");
  BOOST_TEST(matches[0].distance == 2);

  BOOST_TEST(matches[1].name == "down");
  BOOST_TEST(matches[1].distance == 3);

  BOOST_TEST(matches[2].name == "up");
  BOOST_TEST(matches[2].distance == 3);

  BOOST_TEST(matches[3].name == "right");
  BOOST_TEST(matches[3].distance == 5);
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
