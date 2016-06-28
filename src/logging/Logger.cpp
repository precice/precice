#include "Logger.hpp"

#include <boost/log/attributes/mutable_constant.hpp>
#include <boost/log/attributes/named_scope.hpp>
#include <boost/log/attributes/constant.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>

namespace precice {
namespace logging {

Logger::Logger(std::string module)
{
  add_attribute("Module", boost::log::attributes::constant<std::string>(module));
  
  namespace attrs = boost::log::attributes; 
  namespace log = boost::log;
  log::add_common_attributes();
  log::core::get()->add_global_attribute("Scope", attrs::named_scope());
  log::core::get()->add_global_attribute("Rank", attrs::mutable_constant<int>(0));
  log::core::get()->add_global_attribute("Line", attrs::mutable_constant<int>(0));
  log::core::get()->add_global_attribute("File", attrs::mutable_constant<std::string>(""));
  log::core::get()->add_global_attribute("Function", attrs::mutable_constant<std::string>(""));
}

}}
