#pragma once
#define BOOST_LOG_DYN_LINK 1

#include <string>

#include <boost/log/trivial.hpp>
#include <boost/log/sources/severity_logger.hpp>

namespace precice {
namespace logging {

class Logger : public boost::log::sources::severity_logger<boost::log::trivial::severity_level>
{
public:
  explicit Logger(std::string module);
  static void setupLogging();
};
}}// namespace precice, logging
