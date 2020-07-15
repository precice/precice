#include "xml/ValueParser.hpp"
#include <boost/algorithm/string/constants.hpp>
#include <boost/algorithm/string/split.hpp>
#include <ostream>
#include <vector>

namespace precice {
namespace xml {

void readValueSpecific(const std::string &rawValue, Eigen::VectorXd &value)
{
  std::vector<std::string> components;
  boost::split(
      components, rawValue, [](char c) { return c == ';'; }, boost::algorithm::token_compress_on);
  const int size = components.size();
  PRECICE_CHECK(size == 2 || size == 3, "The value \"" << rawValue << "\" is not a 2D or 3D vector.");

  Eigen::VectorXd vec(size);
  for (int i = 0; i != size; ++i) {
    double component{0};
    readValueSpecific(components[i], component);
    vec(i) = component;
  }
  value = vec;
}

} // namespace xml
} // namespace precice
