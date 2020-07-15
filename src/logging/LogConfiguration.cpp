#include "LogConfiguration.hpp"
#include <algorithm>
#include <boost/core/null_deleter.hpp>
#include <boost/log/attributes/mutable_constant.hpp>
#include <boost/log/core.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/support/date_time.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup/console.hpp>
#include <boost/program_options.hpp>
#include <deque>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include "utils/String.hpp"
#include "utils/assertion.hpp"

namespace precice {
namespace logging {

/// A custom formatter that handles the TimeStamp format string
class timestamp_formatter_factory : public boost::log::basic_formatter_factory<char, boost::posix_time::ptime> {
public:
  formatter_type create_formatter(boost::log::attribute_name const &name, args_map const &args)
  {
    namespace expr              = boost::log::expressions;
    args_map::const_iterator it = args.find("format");
    if (it != args.end())
      return expr::stream << expr::format_date_time<boost::posix_time::ptime>(expr::attr<boost::posix_time::ptime>(name), it->second);
    else
      return expr::stream << expr::attr<boost::posix_time::ptime>(name);
  }
};

/// A custom formatter that handles the colorized Severity formatting
class colorized_severity_formatter_factory : public boost::log::formatter_factory<char> {
public:
  formatter_type create_formatter(boost::log::attribute_name const &name, args_map const &args)
  {
    namespace expr = boost::log::expressions;
    auto severity  = expr::attr<boost::log::trivial::severity_level>("Severity");

    return expr::stream
           << expr::if_(severity == boost::log::trivial::severity_level::error)
                  [expr::stream << "\033[31m" //red
                                << "ERROR: "]
           << expr::if_(severity == boost::log::trivial::severity_level::warning)
                  [expr::stream << "\033[36m" //cyan
                                << "WARNING: "]
           << "\033[0m";
  }
};

/// A custom formatter that handles non-colorized Severity formatting
class severity_formatter_factory : public boost::log::formatter_factory<char> {
public:
  formatter_type create_formatter(boost::log::attribute_name const &name, args_map const &args)
  {
    namespace expr = boost::log::expressions;
    auto severity  = expr::attr<boost::log::trivial::severity_level>("Severity");

    return expr::stream
           << expr::if_(severity == boost::log::trivial::severity_level::error)
                  [expr::stream
                   << "ERROR: "]
           << expr::if_(severity == boost::log::trivial::severity_level::warning)
                  [expr::stream
                   << "WARNING: "];
  }
};

/// A simple backends that outputs the message to a stream
/**
 * Rationale: The original text_ostream_backend from boost suffered from the great amount of code that lies
 * between the printing of the message and the endline. This leads to high probability that a process switch 
 * occures and the message is severed from the endline.
 */
class StreamBackend : public boost::log::sinks::text_ostream_backend {
private:
  boost::shared_ptr<std::ostream> _ostream;

public:
  explicit StreamBackend(boost::shared_ptr<std::ostream> ostream)
      : _ostream(ostream) {}

  void consume(boost::log::record_view const &rec, string_type const &formatted_record)
  {
    *_ostream << formatted_record << '\n'
              << std::flush;
  }
};

/// Reads a log file, returns a logging configuration.
LoggingConfiguration readLogConfFile(std::string const &filename)
{
  namespace po = boost::program_options;
  po::options_description desc;
  std::ifstream           ifs(filename);

  po::variables_map vm;

  std::map<std::string, BackendConfiguration> configs;
  try {
    po::parsed_options parsed = parse_config_file(ifs, desc, true);
    po::store(parsed, vm);
    po::notify(vm);
    for (auto const &opt : parsed.options) {
      std::string section = opt.string_key.substr(0, opt.string_key.find('.'));
      std::string key     = opt.string_key.substr(opt.string_key.find('.') + 1);
      configs[section].setOption(key, opt.value[0]);
    }
  } catch (po::error &e) {
    std::cout << "ERROR reading logging configuration: " << e.what() << "\n\n";
    std::exit(-1);
  }

  LoggingConfiguration retVal;

  for (auto const &c : configs)
    if (c.second.enabled)
      retVal.push_back(c.second);

  return retVal;
}

// Default values for filter and format. They are also used from config/LogConfiguration.cpp
const std::string BackendConfiguration::default_filter    = "(%Severity% > debug) and not ((%Severity% = info) and (%Rank% != 0))";
const std::string BackendConfiguration::default_formatter = "(%Rank%) %TimeStamp(format=\"%H:%M:%S\")% [%Module%]:%Line% in %Function%: %ColorizedSeverity%%Message%";
const std::string BackendConfiguration::default_type      = "stream";
const std::string BackendConfiguration::default_output    = "stdout";

void BackendConfiguration::setOption(std::string key, std::string value)
{
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
    enabled = utils::convertStringToBool(value);
  }
}

void setupLogging(LoggingConfiguration configs, bool enabled)
{
  if (_precice_logging_config_lock)
    return;

  namespace bl = boost::log;
  bl::register_formatter_factory("TimeStamp", boost::make_shared<timestamp_formatter_factory>());
  bl::register_formatter_factory("ColorizedSeverity", boost::make_shared<colorized_severity_formatter_factory>());
  bl::register_formatter_factory("Severity", boost::make_shared<severity_formatter_factory>());
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

  for (const auto &config : configs) {
    boost::shared_ptr<StreamBackend> backend;
    if (config.type == "file")
      backend = boost::make_shared<StreamBackend>(boost::shared_ptr<std::ostream>(new std::ofstream(config.output)));
    if (config.type == "stream") {
      if (config.output == "stdout")
        backend = boost::make_shared<StreamBackend>(boost::shared_ptr<std::ostream>(&std::cout, boost::null_deleter()));
      if (config.output == "stderr")
        backend = boost::make_shared<StreamBackend>(boost::shared_ptr<std::ostream>(&std::cerr, boost::null_deleter()));
    }
    PRECICE_ASSERT(backend != nullptr, "The logging backend was not initialized properly. Check your log config.");
    backend->auto_flush(true);
    using sink_t = boost::log::sinks::synchronous_sink<StreamBackend>;
    boost::shared_ptr<sink_t> sink(new sink_t(backend));
    sink->set_formatter(boost::log::parse_formatter(config.format));
    sink->set_filter(boost::log::parse_filter(config.filter));
    boost::log::core::get()->add_sink(sink);
  }
}

void setupLogging(std::string const &logConfigFile)
{
  setupLogging(readLogConfFile(logConfigFile));
}

void setMPIRank(int const rank)
{
  boost::log::attribute_cast<boost::log::attributes::mutable_constant<int>>(boost::log::core::get()->get_global_attributes()["Rank"]).set(rank);
}

void setParticipant(std::string const &participant)
{
  boost::log::attribute_cast<boost::log::attributes::mutable_constant<std::string>>(boost::log::core::get()->get_global_attributes()["Participant"]).set(participant);
}

bool _precice_logging_config_lock{false};

void lockConf()
{
  _precice_logging_config_lock = true;
}

} // namespace logging
} // namespace precice
