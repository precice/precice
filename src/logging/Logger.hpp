#pragma once

#include <string>
#include <boost/log/trivial.hpp>

namespace precice {
namespace logging {

class Logger : public boost::log::sources::severity_logger<boost::log::trivial::severity_level>
{
public:
  explicit Logger(std::string module);
};

}} // namespace precice, logging

// Include LogMacros here, because using it works only together with a Logger
#include "LogMacros.hpp" 

