#include "stacktrace.hpp"

#include <boost/stacktrace.hpp>
#include <sstream>
#include <string>

std::string getStacktrace()
{
  std::ostringstream strm;
  strm << boost::stacktrace::stacktrace();
  return strm.str();
}
