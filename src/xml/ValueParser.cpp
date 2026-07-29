#include "xml/ValueParser.hpp"
#include <boost/algorithm/string/constants.hpp>
#include <boost/algorithm/string/split.hpp>
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

double parseDouble(const std::string &rawValue)
{
  return parseArithmetic<double>(rawValue, "a double");
}
} // namespace

void readValueSpecific(const std::string &rawValue, double &value)
{
  const auto pos = rawValue.find('/');
  if (pos == std::string::npos) {
    value = parseDouble(rawValue);
    return;
  }

  double numerator   = 0.0;
  double denominator = 0.0;
  try {
    numerator   = parseDouble(rawValue.substr(0, pos));
    denominator = parseDouble(rawValue.substr(pos + 1));
  } catch (const std::runtime_error &) {
    // Report the value the user actually wrote, not the offending operand.
    throw std::runtime_error{"Could not fully parse value \"" + rawValue + "\" as a fraction of two doubles."};
  }

  if (denominator == 0.0) {
    throw std::runtime_error{"Could not parse value \"" + rawValue + "\" as a double: the denominator is zero."};
  }
  value = numerator / denominator;
}

void readValueSpecific(const std::string &rawValue, int &value)
{
  value = parseArithmetic<int>(rawValue, "an int");
}

void readValueSpecific(const std::string &rawValue, Eigen::VectorXd &value)
{
  std::vector<std::string> components;
  boost::split(
      components, rawValue, [](char c) { return c == ';'; }, boost::algorithm::token_compress_on);
  const int size = components.size();

  Eigen::VectorXd vec(size);
  for (int i = 0; i != size; ++i) {
    vec(i) = parseDouble(components[i]);
  }
  value = vec;
}

} // namespace precice::xml
