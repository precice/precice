#include "VRML10ParserTest.hpp"
#include "io/impl/VRML10Parser.hpp"
#include "utils/Globals.hpp"
#include "utils/Parallel.hpp"

#ifndef PRECICE_NO_SPIRIT2
#include <boost/spirit/include/qi_parse.hpp>
#include <boost/spirit/include/support_multi_pass.hpp>
namespace spirit = boost::spirit;
#endif // not PRECICE_NO_SPIRIT2

#include <string>
#include <fstream>

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::io::tests::VRML10ParserTest)

namespace precice {
namespace io {
namespace tests {

logging::Logger VRML10ParserTest:: _log ( "precice::io::VRML10ParserTest" );

VRML10ParserTest:: VRML10ParserTest ()
:
  TestCase ( "io::VRML10ParserTest" )
{}

void VRML10ParserTest:: run ()
{
  PRECICE_MASTER_ONLY {
    testMethod ( testParseCube );
  }
}

void VRML10ParserTest:: testParseCube ()
{
# ifndef PRECICE_NO_SPIRIT2
  TRACE();

  for ( int dim=2; dim <= 3; dim++ ){
    std::string file = utils::Globals::getPathToSources() + "/io/tests/ImportVRMLTest";
    DEBUG ( "dim = " << dim );
    if ( dim == 2 ) {
      file += "2D.wrl";
    }
    else {
      assertion ( dim == 3 );
      file += "-Cube.wrl";
    }

    std::ifstream in ( file.c_str() );
    preciceCheck ( in.is_open(), "testParseCube()",
        "Could not open input file: " << file << "!" );

    typedef std::istreambuf_iterator<char> IStreamFileIter;
    typedef spirit::multi_pass<IStreamFileIter> SpiritIter;
    impl::VRML10Parser<SpiritIter> vrmlParser(dim);
    SpiritIter first = spirit::make_default_multi_pass ( IStreamFileIter(in) );
    SpiritIter last = spirit::make_default_multi_pass ( IStreamFileIter() );

    bool success = spirit::qi::phrase_parse ( first, last, vrmlParser, spirit::qi::space );

    if ( first != last ) {
      DEBUG ( "Parsing failed at " << std::string(first, last) );
    }
    validate ( first == last );
    validate ( success );

    if ( dim == 2 ) {
      validateEquals ( vrmlParser.coordinates.size(), 8 );
      validateEquals ( vrmlParser.indices.size(), 6 );
      validateEquals ( vrmlParser.data.size(), 2 );
      validateEquals ( vrmlParser.data[0].name, std::string("Forces") );
      validateEquals ( vrmlParser.data[0].dimensions, 2 );
      validateEquals ( vrmlParser.data[0].values.size(), 8 );
      validateEquals ( vrmlParser.data[1].name, std::string("Velocities") );
      validateEquals ( vrmlParser.data[1].dimensions, 2 );
      validateEquals ( vrmlParser.data[1].values.size(), 8 );
    }
    else {
      assertion ( dim == 3 );
      validateEquals ( vrmlParser.coordinates.size(), 24 );
      validateEquals ( vrmlParser.indices.size(), 36 );
      validateEquals ( vrmlParser.data.size(), 2 );
      validateEquals ( vrmlParser.data[0].name, std::string("Forces") );
      validateEquals ( vrmlParser.data[0].dimensions, 3 );
      validateEquals ( vrmlParser.data[0].values.size(), 24 );
      validateEquals ( vrmlParser.data[1].name, std::string("Velocities") );
      validateEquals ( vrmlParser.data[1].dimensions, 3 );
      validateEquals ( vrmlParser.data[1].values.size(), 24 );
    }
  }
# endif // not PRECICE_NO_SPIRIT2
}

}}} // namespace precice, io, tests
