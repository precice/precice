#include "Logger.hpp"
#include <boost/log/attributes/constant.hpp>
#include <boost/log/attributes/mutable_constant.hpp>
#include <boost/log/attributes/named_scope.hpp>
#include <boost/log/core.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <iosfwd>
#include <utility>

namespace precice {
namespace logging {

/** The implementation of logging::Logger
 *
 * @note The point of using a pimpl for the logger is to remove boost::log from logger.hpp
 */
class Logger::LoggerImpl : public boost::log::sources::severity_logger<boost::log::trivial::severity_level> {
public:
  /** Creates a Boost logger for the said module.
   * @param[in] module the name of the module.
   */
  explicit LoggerImpl(std::string module);
};

Logger::LoggerImpl::LoggerImpl(std::string module)
{
  add_attribute("Module", boost::log::attributes::constant<std::string>(module));

  namespace attrs = boost::log::attributes;
  namespace log   = boost::log;
  log::add_common_attributes();
  log::core::get()->add_global_attribute("Scope", attrs::named_scope());
  log::core::get()->add_global_attribute("Participant", attrs::mutable_constant<std::string>(""));
  log::core::get()->add_global_attribute("Rank", attrs::mutable_constant<int>(0));
  log::core::get()->add_global_attribute("Line", attrs::mutable_constant<int>(0));
  log::core::get()->add_global_attribute("File", attrs::mutable_constant<std::string>(""));
  log::core::get()->add_global_attribute("Function", attrs::mutable_constant<std::string>(""));
}

Logger::Logger(std::string module)
    : _impl(new LoggerImpl{std::move(module)}) {}

// This is required for the std::unique_ptr.
Logger::~Logger() = default;

// Extracts the name saved in the attribute "Module"
Logger::Logger(const Logger &other)
    : Logger{other._impl->get_attributes().find("Module")->second.get_value().extract_or_throw<std::string>()}
{
}

Logger::Logger(Logger &&other) = default;

Logger &Logger::operator=(Logger other)
{
  swap(other);
  return *this;
}

void Logger::swap(Logger &other) noexcept
{
  _impl.swap(other._impl);
}

namespace {
/// Sets the log location in the boost core
void setLogLocation(LogLocation loc)
{
  boost::log::attribute_cast<boost::log::attributes::mutable_constant<int>>(boost::log::core::get()->get_global_attributes()["Line"]).set(loc.line);
  boost::log::attribute_cast<boost::log::attributes::mutable_constant<std::string>>(boost::log::core::get()->get_global_attributes()["File"]).set(loc.file);
  boost::log::attribute_cast<boost::log::attributes::mutable_constant<std::string>>(boost::log::core::get()->get_global_attributes()["Function"]).set(loc.func);
}
} // namespace

void Logger::error(LogLocation loc, const std::string &mess) noexcept
{
  try {
    setLogLocation(loc);
    BOOST_LOG_SEV(*_impl, boost::log::trivial::severity_level::error) << mess;
  } catch (...) {
  }
}

void Logger::warning(LogLocation loc, const std::string &mess) noexcept
{
  try {
    setLogLocation(loc);
    BOOST_LOG_SEV(*_impl, boost::log::trivial::severity_level::warning) << mess;
  } catch (...) {
  }
}

void Logger::info(LogLocation loc, const std::string &mess) noexcept
{
  try {
    setLogLocation(loc);
    BOOST_LOG_SEV(*_impl, boost::log::trivial::severity_level::info) << mess;
  } catch (...) {
  }
}

void Logger::debug(LogLocation loc, const std::string &mess) noexcept
{
  try {
    setLogLocation(loc);
    BOOST_LOG_SEV(*_impl, boost::log::trivial::severity_level::debug) << mess;
  } catch (...) {
  }
}

void Logger::trace(LogLocation loc, const std::string &mess) noexcept
{
  try {
    setLogLocation(loc);
    BOOST_LOG_SEV(*_impl, boost::log::trivial::severity_level::trace) << mess;
  } catch (...) {
  }
}

} // namespace logging
} // namespace precice
