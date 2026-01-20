#include "String.hpp"
#include <Eigen/Dense>
#include <boost/algorithm/string/case_conv.hpp>
#include <memory>
#include <regex>
#include <vector>
#include "boost/algorithm/string/classification.hpp"
#include "boost/algorithm/string/split.hpp"
#include "utils/assertion.hpp"

namespace precice::utils {

std::string wrapText(
    const std::string &text,
    int                linewidth,
    int                indentation)
{
  std::vector<std::string> tokens;
  boost::algorithm::split(tokens, text, boost::algorithm::is_space());
  PRECICE_ASSERT((int) tokens.size() > 0);
  std::string wrapped;
  int         length = 0;
  while (text[length] == ' ') {
    wrapped += ' ';
    length++;
  }
  for (int i = 0; i < static_cast<int>(tokens.size()) - 1; i++) {
    length += static_cast<int>(tokens[i].length());
    wrapped += tokens[i];
    if (length + static_cast<int>(tokens[i + 1].length()) + 1 > linewidth) {
      wrapped += '\n';
      for (int ws = 0; ws < indentation; ws++) {
        wrapped += ' ';
      }
      length = indentation;
    } else {
      wrapped += ' ';
      length += 1;
    }
  }
  wrapped += tokens.back();
  return wrapped;
}

std::string &checkAppendExtension(
    std::string       &filename,
    const std::string &extension)
{
  size_t pos = filename.rfind(extension);
  if ((pos == std::string::npos) || (pos != filename.size() - extension.size())) {
    filename += extension;
  }
  return filename;
}

bool convertStringToBool(std::string const &value)
{
  std::string str = value;
  boost::algorithm::to_lower(str);
  if (str == "1" or str == "yes" or str == "true" or str == "on")
    return true;

  return false;
}

std::string truncate_wstring_to_string(std::wstring wstr, char fill)
{
  std::string converted(wstr.length(), '\0');
  // Buffer for the multibyte representation of a wchar
  std::string mb(MB_CUR_MAX, '\0');
  for (size_t i = 0; i != wstr.length(); ++i) {
    // Converts a wchar to 1-MB_CUR_MAX chars
    std::size_t ret = std::wctomb(&mb[0], wstr[i]);
    converted[i]    = (ret == 1) ? mb.front() : fill;
  }
  return converted;
}

std::size_t editDistance(std::string_view s1, std::string_view s2)
{
  const std::size_t len1 = s1.size(), len2 = s2.size();
  using Matrix = Eigen::Matrix<std::size_t, Eigen::Dynamic, Eigen::Dynamic>;
  Matrix distances(len1 + 1, len2 + 1);

  distances(0, 0) = 0;
  for (std::size_t i = 1; i <= len1; ++i)
    distances(i, 0) = i;
  for (std::size_t i = 1; i <= len2; ++i)
    distances(0, i) = i;

  for (std::size_t i = 1; i <= len1; ++i) {
    for (std::size_t j = 1; j <= len2; ++j) {
      auto deletionCost     = distances(i - 1, j) + 1;
      auto insertionCost    = distances(i, j - 1) + 1;
      auto substitutionCost = distances(i - 1, j - 1) + (s1[i - 1] == s2[j - 1] ? 0 : 1);
      distances(i, j)       = std::min({deletionCost, insertionCost, substitutionCost});
    }
  }

  return distances(len1, len2);
}

bool isKebabStyle(std::string_view sv)
{
  std::regex kebabCaseRegex("^[a-z0-9]+(-[a-z0-9]+)*$");
  return sv.empty() || std::regex_match(sv.begin(), sv.end(), kebabCaseRegex);
}

} // namespace precice::utils
