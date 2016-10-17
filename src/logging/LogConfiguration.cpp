#include "LogConfiguration.hpp"

#include <fstream>
#include <string>

#include <boost/program_options.hpp>

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/attributes/mutable_constant.hpp>
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
  explicit StreamBackend(boost::shared_ptr<std::ostream> ostream) : _ostream(ostream) {}
  
  void consume(boost::log::record_view const& rec, string_type const& formatted_record)
  {
    *_ostream << formatted_record << std::endl << std::flush;
  }
};


/// Reads a log file, returns a logging configuration.
LoggingConfiguration readLogConfFile(std::string filename)
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

  LoggingConfiguration retVal;
    
  for (const auto& c : configs)
    if (c.second.enabled)
      retVal.push_back(c.second);

  return retVal;
}

// Default values for filter and format. They are also used from config/LogConfiguration.cpp
const std::string BackendConfiguration::default_filter = "%Severity% > debug";
const std::string BackendConfiguration::default_formatter = "(%Rank%) %TimeStamp(format=\"%H:%M:%S\")% [%Module%]:%Line% in %Function%: %Message%";
const std::string BackendConfiguration::default_type = "stream";
const std::string BackendConfiguration::default_output = "stdout";

void BackendConfiguration::setOption(std::string key, std::string value)
{
  const std::vector<std::string> trueValues = {"true", "1", "on", "yes"};
  boost::algorithm::to_lower(key);
  if (key == "type") {
    boost::algorithm::to_lower(value);
    type = value;
  }
  if (key == "output")
    output = value;
  if (key == "filter")
    filter = value;
  if (key == "format")
    format = value;
  if (key == "enabled") {
    boost::algorithm::to_lower(value);
    if (std::find(std::begin(trueValues), std::end(trueValues), value) == std::end(trueValues))
      enabled = false;
    else
      enabled = true;
  }
}


void setupLogging(LoggingConfiguration configs, bool enabled)
{
  namespace bl = boost::log;
  bl::register_formatter_factory("TimeStamp", boost::make_shared<timestamp_formatter_factory>());
  bl::register_simple_formatter_factory<bl::trivial::severity_level, char>("Severity");
  bl::register_simple_filter_factory<bl::trivial::severity_level, char>("Severity");

  // Possible, longer output format. Currently unused.
  auto fmtStream =
    bl::expressions::stream
    << "(" << bl::expressions::attr<int>("Rank") << ") "
    << bl::expressions::format_date_time<boost::posix_time::ptime>("TimeStamp", "%H:%M:%S") << " "
    << bl::expressions::attr<std::string>("File") << ":"
    << bl::expressions::attr<int>("Line")
    << " [" << bl::expressions::attr<std::string>("Module") << "] in "
    << bl::expressions::attr<std::string>("Function") << ": "
    << bl::expressions::message;

  // Reset
  bl::core::get()->remove_all_sinks();
  bl::core::get()->reset_filter();

  bl::core::get()->set_logging_enabled(enabled);
  
  // Add the default config
  if (configs.empty())
    configs.emplace_back();
  
  for (const auto& config : configs) {
    boost::shared_ptr<StreamBackend> backend;
    if (config.type == "file")
      backend = boost::make_shared<StreamBackend>(boost::shared_ptr<std::ostream>(new std::ofstream(config.output)));
    if (config.type == "stream") {
      if (config.output == "stdout")
        backend = boost::make_shared<StreamBackend>(boost::shared_ptr<std::ostream>(&std::cout, null_deleter()));
      if (config.output == "stderr")
        backend = boost::make_shared<StreamBackend>(boost::shared_ptr<std::ostream>(&std::cerr, null_deleter()));
    }
    assertion(backend != nullptr, "The logging backend was not initialized properly. Check your log config.");
    backend->auto_flush(true);
    using sink_t =  boost::log::sinks::synchronous_sink<StreamBackend>;          
    boost::shared_ptr<sink_t> sink(new sink_t(backend));
    sink->set_formatter(boost::log::parse_formatter(config.format));
    sink->set_filter(boost::log::parse_filter(config.filter));
    boost::log::core::get()->add_sink(sink);
  }    
}


void setupLogging(std::string logConfigFile)
{
  setupLogging(readLogConfFile(logConfigFile));
}


void setMPIRank(const int rank) {
  boost::log::attribute_cast<boost::log::attributes::mutable_constant<int>>(boost::log::core::get()->get_global_attributes()["Rank"]).set(rank);
}

}} // namespace precice, logging
