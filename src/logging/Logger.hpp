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
  
};

void setupLogging(std::string logConfigFile = "log.conf");

/// Sets the current MPI rank as a logging attribute
void setMPIRank(const int rank);

}}// namespace precice, logging
