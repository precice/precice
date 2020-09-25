#include "xml/ValueParser.hpp"
#include <boost/algorithm/string/constants.hpp>
#include <boost/algorithm/string/split.hpp>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace precice {
namespace xml {

namespace {
constexpr static const char *PARSING_LOCALE = "en_US.UTF-8";

double parseDouble(const std::string &rawValue)
{
  std::istringstream iss{rawValue};
  try {
    iss.imbue(std::locale(PARSING_LOCALE));
  } catch (...) {
  }
  double value;
  iss >> value;
  if (!iss.eof()) {
    throw std::runtime_error{"Could not fully parse value \"" + rawValue + "\" as a double."};
  }
  return value;
}
} // namespace

void readValueSpecific(const std::string &rawValue, double &value)
{
  if (rawValue.find('/') != std::string::npos) {
    std::string left  = rawValue.substr(0, rawValue.find('/'));
    std::string right = rawValue.substr(rawValue.find('/') + 1, rawValue.size() - rawValue.find('/') - 1);

    value = parseDouble(left) / parseDouble(right);
  } else {
    value = parseDouble(rawValue);
  }
}

void readValueSpecific(const std::string &rawValue, int &value)
{
  std::istringstream iss{rawValue};
  try {
    iss.imbue(std::locale(PARSING_LOCALE));
  } catch (...) {
  }
  iss >> value;
  if (!iss.eof()) {
    throw std::runtime_error{"Could not fully parse value \"" + rawValue + "\" as an int."};
  }
}

void readValueSpecific(const std::string &rawValue, Eigen::VectorXd &value)
{
  std::vector<std::string> components;
  boost::split(
      components, rawValue, [](char c) { return c == ';'; }, boost::algorithm::token_compress_on);
  const int size = components.size();
  if (size < 2 || size > 3) {
    throw std::runtime_error{"The value \"" + rawValue + "\" is not a 2D or 3D vector."};
  }

  Eigen::VectorXd vec(size);
  for (int i = 0; i != size; ++i) {
    vec(i) = parseDouble(components[i]);
  }
  value = vec;
}

} // namespace xml
} // namespace precice
