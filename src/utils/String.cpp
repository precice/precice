#include "String.hpp"
#include <boost/algorithm/string/case_conv.hpp>
#include <memory>
#include <vector>
#include "boost/algorithm/string/classification.hpp"
#include "boost/algorithm/string/split.hpp"
#include "utils/assertion.hpp"

namespace precice {
namespace utils {

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
  for (int i = 0; i < (int) tokens.size() - 1; i++) {
    length += (int) tokens[i].length();
    wrapped += tokens[i];
    if (length + (int) tokens[i + 1].length() + 1 > linewidth) {
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
    std::string &      filename,
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

} // namespace utils
} // namespace precice
