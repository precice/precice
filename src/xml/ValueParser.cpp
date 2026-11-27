#include "xml/ValueParser.hpp"
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace precice::xml {

namespace {
constexpr static const char *PARSING_LOCALE = "en_US.UTF-8";

/// Extracts a single arithmetic value, requiring the entire string to be consumed.
template <typename T>
T parseArithmetic(const std::string &rawValue, const char *typeName)
{
  std::istringstream iss{rawValue};
  try {
    iss.imbue(std::locale(PARSING_LOCALE));
  } catch (...) {
  }
  T value{};
  iss >> value;
  // An empty or whitespace-only string makes the sentry fail, which sets eofbit alongside
  // failbit and leaves value untouched. Checking eof() alone therefore misses these.
  if (iss.fail() || !iss.eof()) {
    throw std::runtime_error{"Could not fully parse value \"" + rawValue + "\" as " + typeName + "."};
  }
  return value;
}
} // namespace

void readValueSpecific(const std::string &rawValue, double &value)
{
  value = parseArithmetic<double>(rawValue, "a double");
}

void readValueSpecific(const std::string &rawValue, int &value)
{
  value = parseArithmetic<int>(rawValue, "an int");
}

void readValueSpecific(const std::string &rawValue, Eigen::VectorXd &value)
{
  std::vector<std::string> components;
  boost::split(
      components, rawValue, [](char c) { return c == ';'; });
  const Eigen::Index size = static_cast<Eigen::Index>(components.size());

  Eigen::VectorXd vec(size);
  for (Eigen::Index i = 0; i != size; ++i) {
    boost::trim(components[i]);
    vec(i) = parseArithmetic<double>(components[i], "a double");
  }
  value = vec;
}

} // namespace precice::xml
