#include "Logger.hpp"

#include <iostream>
#include <fstream>
#include <string>
using namespace std;

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

// #include <boost/log/utility/manipulators/add_value.hpp>

#include <boost/log/attributes/named_scope.hpp>
#include <boost/log/attributes/mutable_constant.hpp>

//#include <boost/log/sources/severity_logger.hpp>

#include <boost/log/expressions/keyword.hpp>
#include <boost/log/expressions/formatters/named_scope.hpp>

#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/log/utility/setup/console.hpp>
#include <boost/log/utility/setup/from_stream.hpp>

#include <boost/log/support/date_time.hpp>


namespace precice {
namespace logging {

Logger::Logger(std::string module)
{
  add_attribute("Module", boost::log::attributes::constant<std::string>(module));
}

void setupLogging()
{
  using namespace boost::log;
  add_common_attributes();
  core::get()->add_global_attribute("Scope", attributes::named_scope());
  core::get()->add_global_attribute("Rank", attributes::mutable_constant<int>(0));
  core::get()->add_global_attribute("Line", attributes::mutable_constant<int>(0));
  core::get()->add_global_attribute("File", attributes::mutable_constant<std::string>(""));
  core::get()->add_global_attribute("Function", attributes::mutable_constant<std::string>(""));

  register_simple_formatter_factory<trivial::severity_level, char>("Severity");
  register_simple_filter_factory<trivial::severity_level, char>("Severity");

  std::string format = "%LineID%:Code: [%TimeStamp%] [Rank: "
    "%Rank%][%Scope%:%Line%] [%Severity%]: %Message%";

  auto fmtStream =
    expressions::stream
    << "(" 
    << expressions::attr<int>("Rank")
    << ") "
    << expressions::format_date_time<boost::posix_time::ptime>("TimeStamp", "%H:%M:%S")
    << " "
    << expressions::attr<std::string>("File")
    << ":"
    << expressions::attr<int>("Line")
    << " ["
    << expressions::attr<std::string>("Module")
    << "] in "
    << expressions::attr<std::string>("Function") 
    << ": "
    << expressions::message << std::endl;
    //<< "LineID: " << expressions::attr<unsigned int>("LineID") << std::endl
    //<< ", ThreadID: " << expressions::attr<attributes::current_thread_id::value_type>("ThreadID") << " "
    //<< ", ProcessID: " << expressions::attr<attributes::current_process_id::value_type>("ProcessID") << std::endl
    //<< ", Severity: " << trivial::severity
    //<< "Scope: " << expressions::format_named_scope("Scope", keywords::format = "%n in %f:%l)") << std::endl

  boost::log::add_file_log("sample.log", keywords::format = fmtStream); // state namespace here for
                                                                        // consisitency with
                                                                        // console_log
  boost::log::add_console_log(std::cout, keywords::format = fmtStream); // explicitly state namespace
                                                                        // here to resolve some
                                                                        // ambiguity issue

  std::ifstream file("log.conf");
  // settings conf = parse_settings(file);
  // cout << conf["Core"]["Filter"] << endl;

  if (file.is_open())
    init_from_stream(file);
}

void setMPIRank(const int rank){
  boost::log::attribute_cast<boost::log::attributes::mutable_constant<int>>(boost::log::core::get()->get_global_attributes()["Rank"]).set(rank);
}
}}// namespace precice, logging
