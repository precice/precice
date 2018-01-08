#ifndef PRECICE_NO_SPIRIT2

#include <boost/spirit/include/qi_parse.hpp>
#include <boost/spirit/include/support_multi_pass.hpp>
#include <fstream>
#include <string>
#include "io/impl/VRML10Parser.hpp"
#include "testing/Testing.hpp"
#include "utils/Globals.hpp"

BOOST_AUTO_TEST_SUITE(VRML10ParserTests)

using namespace precice;
using namespace precice::io;
namespace spirit = boost::spirit;

BOOST_AUTO_TEST_CASE(ParseCube)
{
  for (int dim = 2; dim <= 3; dim++) {
    std::string file = utils::getPathToSources() + "/io/tests/ImportVRMLTest";
    if (dim == 2) {
      file += "2D.wrl";
    } else {
      assertion(dim == 3);
      file += "-Cube.wrl";
    }

    std::ifstream in(file.c_str());
    BOOST_REQUIRE_MESSAGE(in.is_open(), "Could not open input file: " << file << "!");

    using IStreamFileIter = std::istreambuf_iterator<char>;
    using SpiritIter      = spirit::multi_pass<IStreamFileIter>;
    impl::VRML10Parser<SpiritIter> vrmlParser(dim);
    SpiritIter                     first = spirit::make_default_multi_pass(IStreamFileIter(in));
    SpiritIter                     last  = spirit::make_default_multi_pass(IStreamFileIter());

    bool success = spirit::qi::phrase_parse(first, last, vrmlParser, spirit::qi::space);

    BOOST_CHECK_MESSAGE(first == last, "Parsing failed at " << std::string(first, last));
    BOOST_TEST(success);

    if (dim == 2) {
      BOOST_TEST(vrmlParser.coordinates.size() == 8);
      BOOST_TEST(vrmlParser.indices.size() == 6);
      BOOST_TEST(vrmlParser.data.size() == 2);
      BOOST_TEST(vrmlParser.data[0].name == "Forces");
      BOOST_TEST(vrmlParser.data[0].dimensions == 2);
      BOOST_TEST(vrmlParser.data[0].values.size() == 8);
      BOOST_TEST(vrmlParser.data[1].name == "Velocities");
      BOOST_TEST(vrmlParser.data[1].dimensions == 2);
      BOOST_TEST(vrmlParser.data[1].values.size() == 8);
    } else {
      assertion(dim == 3);
      BOOST_TEST(vrmlParser.coordinates.size() == 24);
      BOOST_TEST(vrmlParser.indices.size() == 36);
      BOOST_TEST(vrmlParser.data.size() == 2);
      BOOST_TEST(vrmlParser.data[0].name == "Forces");
      BOOST_TEST(vrmlParser.data[0].dimensions == 3);
      BOOST_TEST(vrmlParser.data[0].values.size() == 24);
      BOOST_TEST(vrmlParser.data[1].name == "Velocities");
      BOOST_TEST(vrmlParser.data[1].dimensions == 3);
      BOOST_TEST(vrmlParser.data[1].values.size() == 24);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END() // IOTests

#endif // not PRECICE_NO_SPIRIT2
