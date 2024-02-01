#include <boost/log/attributes/constant.hpp>
#include <boost/log/attributes/function.hpp>
#include <boost/log/attributes/named_scope.hpp>
#include <boost/log/attributes/timer.hpp>
#include <boost/log/sources/severity_feature.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <utility>
#include <utils/assertion.hpp>

#include "logging/LogConfiguration.hpp"
#include "logging/Logger.hpp"

namespace precice::logging {

/// A feature that adds \ref LogLocation info to log attributes
template <class BaseT>
class precice_feature : public BaseT {
public:
  using char_type       = typename BaseT::char_type;
  using threading_model = typename BaseT::threading_model;

  using open_record_lock = typename boost::log::strictest_lock<
      boost::lock_guard<threading_model>,
      typename BaseT::open_record_lock,
      typename BaseT::add_attribute_lock,
      typename BaseT::remove_attribute_lock>::type;

protected:
  template <typename ArgsT>
  boost::log::record open_record_unlocked(ArgsT const &args);
};

/// Provides keywords to pass \ref LogLocation related info to the precice_feature
/// @see PRECICE_LOG_IMPL
namespace keywords {
BOOST_PARAMETER_KEYWORD(line_ns, line);
BOOST_PARAMETER_KEYWORD(file_ns, file);
BOOST_PARAMETER_KEYWORD(func_ns, func);
} // namespace keywords

/// Adds log attributes to the current logger based on the \ref LogLocation info
template <typename BaseT>
template <typename ArgsT>
boost::log::record precice_feature<BaseT>::open_record_unlocked(ArgsT const &args)
{
  // Extract the named argument from the parameters pack
  int         line_value = args[keywords::line | -1];
  std::string file_value = args[keywords::file | std::string()];
  std::string func_value = args[keywords::func | std::string()];

  PRECICE_ASSERT(line_value >= 0);
  PRECICE_ASSERT(!file_value.empty());
  PRECICE_ASSERT(!func_value.empty());

  // Process extracted tags
  using boost::log::attributes::constant;
  boost::log::attribute_set &attrs = BaseT::attributes();
  {
    [[maybe_unused]] auto res = BaseT::add_attribute_unlocked("Line", constant<int>(line_value));
    PRECICE_ASSERT(res.second);
  }
  {
    [[maybe_unused]] auto res = BaseT::add_attribute_unlocked("File", constant<std::string>(std::move(file_value)));
    PRECICE_ASSERT(res.second);
  }
  {
    [[maybe_unused]] auto res = BaseT::add_attribute_unlocked("Function", constant<std::string>(std::move(func_value)));
    PRECICE_ASSERT(res.second);
  }

  // Forward the call to the base feature
  auto &&record = BaseT::open_record_unlocked(args);

  // cleanup
  attrs.erase("Line");
  attrs.erase("File");
  attrs.erase("Function");

  return std::move(record);
}

/// Required by \ref BoostLogger to use \ref precice_feature
struct precice_log : public boost::mpl::quote1<precice_feature> {
};

/// The boost logger that combines required featrues
template <class BaseLogger>
using BoostLogger = boost::log::sources::basic_composite_logger<
    char,
    BaseLogger,
    boost::log::sources::single_thread_model,
    boost::log::sources::features<
        boost::log::sources::severity<boost::log::trivial::severity_level>,
        precice_log>>;

/** The implementation of logging::Logger
 *
 * @note The point of using a pimpl for the logger is to remove boost::log from logger.hpp
 */
class Logger::LoggerImpl : public BoostLogger<Logger::LoggerImpl> {
public:
  /** Creates a Boost logger for the said module.
   * @param[in] module the name of the module.
   */
  explicit LoggerImpl(std::string_view module);
};

/// Registers attributes that don't depend on the \ref LogLocation
Logger::LoggerImpl::LoggerImpl(std::string_view module)
{
  namespace attrs = boost::log::attributes;

  // The preCICE attribute marks attribute values that originate from preCICE
  add_attribute("preCICE", attrs::constant<bool>(true));
  add_attribute("Participant", attrs::make_function<std::string>([] { return getGlobalLoggingConfig().participant; }));
  add_attribute("Rank", attrs::make_function<int>([] { return getGlobalLoggingConfig().rank; }));

  // Constant attributes
  add_attribute("Module", attrs::constant<std::string>(std::string{module}));

  // Ensure, boost-provided attributes are present
  add_attribute("TimeStamp", attrs::local_clock());
  add_attribute("Runtime", attrs::timer());
  add_attribute("Scope", attrs::named_scope());
}

Logger::Logger(std::string_view module)
    : _impl(new LoggerImpl{module}) {}

/// This is required for the std::unique_ptr.
Logger::~Logger() = default;

/// Implements a deep copy of the implementation
Logger::Logger(const Logger &other)
    : _impl(std::make_unique<LoggerImpl>(*other._impl))
{
}

Logger::Logger(Logger &&other) noexcept = default;

Logger &Logger::operator=(Logger other)
{
  swap(other);
  return *this;
}

void Logger::swap(Logger &other) noexcept
{
  _impl.swap(other._impl);
}

/// Convenience macro that automatically passes \ref LogLocation info to the logger
#define PRECICE_LOG_IMPL(lg, sev, loc) \
  BOOST_LOG_STREAM_WITH_PARAMS((lg), (boost::log::keywords::severity = (sev))(keywords::line = (loc.line))(keywords::file = (loc.file))(keywords::func = (loc.func)))

void Logger::error(LogLocation loc, std::string_view mess) noexcept
{
  try {
    PRECICE_LOG_IMPL(*_impl, boost::log::trivial::severity_level::error, loc) << mess;
  } catch (...) {
  }
}

void Logger::warning(LogLocation loc, std::string_view mess) noexcept
{
  try {
    PRECICE_LOG_IMPL(*_impl, boost::log::trivial::severity_level::warning, loc) << mess;
  } catch (...) {
  }
}

void Logger::info(LogLocation loc, std::string_view mess) noexcept
{
  try {
    PRECICE_LOG_IMPL(*_impl, boost::log::trivial::severity_level::info, loc) << mess;
  } catch (...) {
  }
}

void Logger::debug(LogLocation loc, std::string_view mess) noexcept
{
  try {
    PRECICE_LOG_IMPL(*_impl, boost::log::trivial::severity_level::debug, loc) << mess;
  } catch (...) {
  }
}

void Logger::trace(LogLocation loc, std::string_view mess) noexcept
{
  try {
    PRECICE_LOG_IMPL(*_impl, boost::log::trivial::severity_level::trace, loc) << mess;
  } catch (...) {
  }
}

#undef PRECICE_LOG_IMPL

} // namespace precice::logging
