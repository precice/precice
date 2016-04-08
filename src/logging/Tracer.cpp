#include "Tracer.hpp"
#include <boost/log/attributes/mutable_constant.hpp>

namespace precice {
namespace logging {

Tracer::Tracer
(
  Logger &log,
  std::string function,
  std::string file,
  long line
  )
  :
  _log(log),
  _function(function),
  _file(file),
  _line(line)
{}

Tracer::~Tracer()
{
  using namespace boost::log;

  attribute_cast<attributes::mutable_constant<int>>(core::get()->get_global_attributes()["Line"]).set(_line);

  attribute_cast<attributes::mutable_constant<std::string>>(core::get()->get_global_attributes()["File"]).set(_file);

  attribute_cast<attributes::mutable_constant<std::string>>(core::get()->get_global_attributes()["Function"]).set(_function);

  BOOST_LOG_SEV(_log, trivial::severity_level::trace) << "Leaving " << _function;
}

}} // namespace precice,logging
