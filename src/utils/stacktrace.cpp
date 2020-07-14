#include "stacktrace.hpp"
#include <boost/stacktrace.hpp>
#include <exception>
#include <sstream>
#include <string>

std::string getStacktrace()
{
  std::ostringstream strm;
  try {
    strm << boost::stacktrace::stacktrace();
  } catch (const std::exception &e) {
    strm << "Stacktrace failed: " << e.what();
  }
  return strm.str();
}
