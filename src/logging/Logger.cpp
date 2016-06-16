#include "Logger.hpp"

#include <fstream>
#include <string>

#include <boost/program_options.hpp>

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/attributes/named_scope.hpp>
#include <boost/log/attributes/mutable_constant.hpp>
#include <boost/log/expressions/keyword.hpp>
#include <boost/log/expressions/formatters/named_scope.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/log/utility/setup/console.hpp>
#include <boost/log/support/date_time.hpp>

#include "utils/assertion.hpp"

namespace precice {
namespace logging {

/// Function object that does nothing
/* boost 1.55 ships the empty_deleter that was later renamend to null_deleter. boost 1.54 includes neither.
 * This functor is used together with boost::shared_ptr to avoid deleting std::cout/cerr.
 * It can be removed and replaced by empty_deleter if boost 1.54 does not need to be supported anymore.
 * Omitting the deleter alltogether causes a double-free or corruption.
 */
struct null_deleter
{
  template<typename T> void operator() (T*) const {}
};



/// A custom formatter that handle the TimeStamp format string
class timestamp_formatter_factory :
    public boost::log::basic_formatter_factory<char, boost::posix_time::ptime>
{
public:
  formatter_type create_formatter(boost::log::attribute_name const& name, args_map const& args)
  {
    namespace expr = boost::log::expressions;
    args_map::const_iterator it = args.find("format");
    if (it != args.end())
      return expr::stream <<
        expr::format_date_time<boost::posix_time::ptime>(expr::attr<boost::posix_time::ptime>(name), it->second);
    else
      return expr::stream <<
        expr::attr<boost::posix_time::ptime>(name);
  }
};

/// Holds the configuration for one logging backend and takes care of default values.
struct BackendConfiguration
{
  std::string type;
  std::string output = "stdout";
  boost::log::filter filter = boost::log::parse_filter("%Severity% > debug");
  boost::log::basic_formatter<char> format = boost::log::parse_formatter("(%Rank%) %TimeStamp(format=\"%H:%M:%S\")% [%Module%]:%LINE% in %Function%: %Message%");

  void setOption(std::string key, std::string value)
  {
    boost::algorithm::to_lower(key);
    if (key == "type") {
      boost::algorithm::to_lower(value);
      type = value;
    }
    if (key == "output")
      output = value;
    if (key == "filter")
      filter = boost::log::parse_filter(value);
    if (key == "format")
      format = boost::log::parse_formatter(value);
  }    
};


/// A simple backends that outputs the message to a stream
/* Rationale: The original test_ostream_backend from boost suffered from the great amount of code that lies 
 * between the printing of the message and the endline. This leads to high probability that a process switch 
 * occures and the message is severed from the endline.
 */
class StreamBackend :
    public boost::log::sinks::text_ostream_backend
{
private:
  boost::shared_ptr<std::ostream> _ostream;
  
public:
  StreamBackend(boost::shared_ptr<std::ostream> ostream) : _ostream(ostream) {}
  
  void consume(boost::log::record_view const& rec, string_type const& formatted_record)
  {
    *_ostream << formatted_record << std::endl << std::flush;
  }
};


Logger::Logger(std::string module)
{
  add_attribute("Module", boost::log::attributes::constant<std::string>(module));
}


std::map<std::string, BackendConfiguration> readLogConfFile(std::string filename)
{
  namespace po = boost::program_options;
  po::options_description desc;
  std::ifstream ifs(filename);
    
  po::variables_map vm;

  std::map<std::string, BackendConfiguration> configs;
  try {
    po::parsed_options parsed = parse_config_file(ifs, desc, true);
    po::store(parsed, vm);
    po::notify(vm);
    for (const auto& opt : parsed.options) {
      std::string section = opt.string_key.substr(0, opt.string_key.find("."));
      std::string key = opt.string_key.substr(opt.string_key.find(".")+1);
      configs[section].setOption(key, opt.value[0]);        
    }
  }
  catch (po::error& e) {
    std::cout << "ERROR reading logging configuration: " << e.what() << "\n\n";
    std::exit(-1);
  }
  return configs;  
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

  // Possible, longer output format. Currently unused.
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
     << expressions::message;
  
  auto configs = readLogConfFile(logConfigFile);

  if (configs.empty()) {
    configs["DefaultBackend"].type = "stream";
    configs["DefaultBackend"].output = "stdout";
  }
  for (const auto& config : configs) {
    std::cout << "Initialized Logging Backend: " << config.first << std::endl; // use std::cout here, because logging isn't initialized yet
    boost::shared_ptr<StreamBackend> backend;
    if (config.second.type == "file")
      backend = boost::make_shared<StreamBackend>(boost::shared_ptr<std::ostream>(new std::ofstream(config.second.output)));
    if (config.second.type == "stream") {
      if (config.second.output == "stdout")
        backend = boost::make_shared<StreamBackend>(boost::shared_ptr<std::ostream>(&std::cout, null_deleter()));
      if (config.second.output == "stderr")
        backend = boost::make_shared<StreamBackend>(boost::shared_ptr<std::ostream>(&std::cerr, null_deleter()));
    }
    assertion(backend != nullptr, "The logging backend was not initialized properly. Check your " + logConfigFile);
    backend->auto_flush(true);
    using sink_t =  boost::log::sinks::synchronous_sink<StreamBackend>;          
    boost::shared_ptr<sink_t> sink(new sink_t(backend));
    sink->set_formatter(config.second.format);
    sink->set_filter(config.second.filter);
    boost::log::core::get()->add_sink(sink);
  }    
}

void setMPIRank(const int rank) {
  boost::log::attribute_cast<boost::log::attributes::mutable_constant<int>>(boost::log::core::get()->get_global_attributes()["Rank"]).set(rank);
}

}} // namespace precice, logging
