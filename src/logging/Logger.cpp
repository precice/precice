#include "Logger.hpp"

#include <fstream>
#include <string>

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
#include <boost/log/utility/setup/from_settings.hpp>
#include <boost/log/utility/setup/from_stream.hpp>
#include <boost/log/utility/setup/settings.hpp>

#include <boost/log/support/date_time.hpp>


namespace precice {
namespace logging {

class timestamp_formatter_factory :
    public boost::log::basic_formatter_factory<char, boost::posix_time::ptime>
{
public:
  formatter_type create_formatter(boost::log::attribute_name const& name, args_map const& args)
  {
    args_map::const_iterator it = args.find("format");
    if (it != args.end())
      return boost::log::expressions::stream <<
        boost::log::expressions::format_date_time<boost::posix_time::ptime>(boost::log::expressions::attr<boost::posix_time::ptime>(name), it->second);
    else
      return boost::log::expressions::stream <<
        boost::log::expressions::attr<boost::posix_time::ptime>(name);
  }
};

Logger::Logger(std::string module)
{
  add_attribute("Module", boost::log::attributes::constant<std::string>(module));
}

void setupLogging(std::string logConfigFile)
{
  using namespace boost::log;
  register_formatter_factory("TimeStamp", boost::make_shared<timestamp_formatter_factory>());
  add_common_attributes();
  core::get()->add_global_attribute("Scope", attributes::named_scope());
  core::get()->add_global_attribute("Rank", attributes::mutable_constant<int>(0));
  core::get()->add_global_attribute("Line", attributes::mutable_constant<int>(0));
  core::get()->add_global_attribute("File", attributes::mutable_constant<std::string>(""));
  core::get()->add_global_attribute("Function", attributes::mutable_constant<std::string>(""));

  register_simple_formatter_factory<trivial::severity_level, char>("Severity");
  register_simple_filter_factory<trivial::severity_level, char>("Severity");

  std::string format = "(%Rank%) %TimeStamp(format=\"%H:%M:%S\")% %File%:%Line%[%Module%] in %Function%: %Message%";
  
  settings setts;

  setts["Core.Filter"] = "%Severity% > debug";
  setts["Core.DisableLogging"] = false;

  // Subsections can be referred to with a single path
  setts["Sinks.Console.Destination"] = "Console";
  setts["Sinks.Console.Filter"] = "%Severity% > debug";
  setts["Sinks.Console.AutoFlush"] = true;
  setts["Sinks.Console.Format"] = format;

  // ...as well as the individual parameters
/*  setts["Sinks.File.Destination"] = "TextFile";
  setts["Sinks.File.FileName"] = "sample.log";
  setts["Sinks.File.AutoFlush"] = true;
  setts["Sinks.File.RotationSize"] = 10 * 1024 * 1024; // 10 MiB
  setts["Sinks.File.Format"] = format;
*/
  
//alternative setting of log format
  /*
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
     << expressions::message; //<< std::endl;
  */
//Additional possibilities for debugging output
//expressions::attr<unsigned int>("LineID")
//expressions::attr<attributes::current_thread_id::value_type>("ThreadID")
//expressions::attr<attributes::current_process_id::value_type>("ProcessID")
//trivial::severity
//expressions::format_named_scope("Scope", keywords::format = "%n in %f:%l)")

  
  //boost::log::add_file_log("sample.log", keywords::format = fmtStream); // state namespace here for
  // consisitency with
  // console_log
  //boost::log::add_console_log(std::cout, keywords::format = fmtStream); // explicitly state namespace
  // here to resolve some
  // ambiguity issue

  //if config file exists only entries in the file overrides our standard config only
  std::ifstream file(logConfigFile);
  // settings conf = parse_settings(file);
  // cout << conf["Core"]["Filter"] << endl;

  if (file.is_open()){
    init_from_stream(file);
  } else {
    init_from_settings(setts);
  }
}

void setMPIRank(const int rank){
  boost::log::attribute_cast<boost::log::attributes::mutable_constant<int>>(boost::log::core::get()->get_global_attributes()["Rank"]).set(rank);
}
}}// namespace precice, logging
