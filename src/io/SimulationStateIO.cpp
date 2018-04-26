#include "SimulationStateIO.hpp"
#include "utils/Globals.hpp"
#ifndef PRECICE_NO_SPIRIT2
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/qi_parse.hpp>
#include <boost/spirit/include/support_multi_pass.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
namespace spirit = boost::spirit;
namespace qi = boost::spirit::qi;
namespace phoenix = boost::phoenix;
#endif // not PRECICE_NO_SPIRIT2
#include <fstream>

namespace precice {
namespace io {

//const std::string& SimulationStateIO:: standardFileName()
//{
//  static std::string name = "precice_";
//  return name;
//}

const std::string& SimulationStateIO:: standardFileExtension()
{
  static std::string ext = ".txt";
  return ext;
}

SimulationStateIO:: SimulationStateIO
(
  const std::string& file)
:
  _file(file)
{
  assertion(file != std::string(""));
}

void SimulationStateIO:: writeState
(
  double   time,
  int      timestep,
  long int numberAdvanceCalls)
{
  std::ofstream out;
  out.open(_file.c_str());
  if (not out){
    ERROR("Could not open file \"" << _file
                 << "\" for writing the simulation state!");
  }
  out << "preCICE simulation state" << std::endl
      << std::endl
      << "time " << time << std::endl
      << "timestep " << timestep << std::endl
      << "advancecalls " << numberAdvanceCalls << std::endl;
  out.close();
}

void SimulationStateIO:: readState
(
  double&   time,
  int&      timestep,
  long int& numberAdvanceCalls)
{
# ifndef PRECICE_NO_SPIRIT2
  std::ifstream in(_file.c_str());
  if (not in.is_open()){
    ERROR("Could not open file " << _file
                 << " to read simulation state!");
  }
  // Wrap input file stream into a multi pass iterator (requirement)
  typedef std::istreambuf_iterator<char> FileIter;
  typedef spirit::multi_pass<FileIter> MultiPassFileIter;
  MultiPassFileIter first = spirit::make_default_multi_pass(FileIter(in));
  MultiPassFileIter last = spirit::make_default_multi_pass(FileIter());
  // Parse VRML file and check validity
  qi::rule<MultiPassFileIter, qi::space_type> parser =
     qi::lit("preCICE simulation state")
     >> qi::lit("time") >> qi::double_[phoenix::ref(time) = qi::_1]
     >> qi::lit("timestep") >> qi::int_[phoenix::ref(timestep) = qi::_1]
     >> qi::lit("advancecalls") >> qi::long_[phoenix::ref(numberAdvanceCalls) = qi::_1];
  bool success = qi::phrase_parse(first, last, parser, qi::space);
  if((!success) || (first != last)) {
     ERROR(
                  "Parsing of file " << _file << " failed! Left over: "
                  << std::endl << std::string(first, last));
  }
# else // not PRECICE_NO_SPIRIT2
  ERROR("Needs Boost.Spirit V2.0!");
# endif // PRECICE_NO_SPIRIT2
}

}} // namespace precice, io
