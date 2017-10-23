#include <string>
#include "testing/Testing.hpp"
#include "utils/String.hpp"

using namespace precice::utils;

BOOST_AUTO_TEST_SUITE(UtilsTests)

BOOST_AUTO_TEST_CASE(StringWrap)
{
  std::string text("123456 1234567 12345678");
  std::string wrapped = wrapText(text, 4, 0);
  BOOST_TEST(wrapped == std::string("123456\n1234567\n12345678"));

  text = "1234 1234 1234";
  wrapped = wrapText(text, 5, 0);
  BOOST_TEST(wrapped == std::string("1234\n1234\n1234"));

  text = "1234 123 5";
  wrapped = wrapText(text, 5, 0);
  BOOST_TEST(wrapped == std::string("1234\n123 5"));

  text = "1234 1234 1234";
  wrapped = wrapText(text, 7, 2);
  BOOST_TEST(wrapped == std::string("1234\n  1234\n  1234"));

  text = "1234 1234 1234";
  wrapped = wrapText(text, 5, 2);
  BOOST_TEST(wrapped == std::string("1234\n  1234\n  1234"));

  text = "12345678 1234 1 1234";
  wrapped = wrapText(text, 8, 2);
  BOOST_TEST(wrapped == std::string("12345678\n  1234 1\n  1234"));
}

BOOST_AUTO_TEST_CASE(StringAppendExtension)
{
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

BOOST_AUTO_TEST_SUITE_END()
