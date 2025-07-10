#pragma once

#include <algorithm>
#include <array>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

namespace precice::utils {

/// Utility class to build a string from C functions with output pointers and static maximum length
template <int MAX>
class StringMaker {
public:
  StringMaker()
  {
    clear();
  }

  void clear()
  {
    _data.fill('\0');
  }

  char *data()
  {
    return _data.data();
  }

  /** constructs a string from the buffer
   *
   *  The returned string ends at the fill NULL char.
   */
  std::string str() const
  {
    return std::string(_data.data());
  }

private:
  std::array<char, MAX + 1> _data;
};

std::string wrapText(
    const std::string &text,
    int                linewidth,
    int                indentation);

/**
 * @brief Checks if filename has the given extension, if not appends it.
 *
 * @return filename with extension.
 */
std::string &checkAppendExtension(std::string &filename, const std::string &extension);

/// Evaluates a string to find out if it represents a bool.
/**
 * Returns True if string is yes, true, 1 or on. Otherwise False.
 * This function is case-insensitive.
 */
bool convertStringToBool(std::string const &value);

/** Converts a wstring to a string by truncating non-extended ascii characters.
 *
 * @param[in] wstr the wide string to convert
 * @param[in] fill the fill char to replace multibyte characters with
 * @return the converted string
 */
std::string truncate_wstring_to_string(std::wstring wstr, char fill = '#');

/// Turns stream-like code into a std::string.
#define PRECICE_AS_STRING(message) [&] { \
  std::ostringstream oss;                \
  oss << message;                        \
  return oss.str();                      \
}()

std::size_t editDistance(std::string_view s1, std::string_view s2);

struct StringMatch {
  std::string name;
  std::size_t distance;

  bool operator<(const StringMatch &other) const
  {
    if (distance == other.distance) {
      return name < other.name;
    } else {
      return distance < other.distance;
    }
  }
};

template <class Container>
std::vector<StringMatch> computeMatches(std::string_view given, const Container &expected)
{
  std::vector<StringMatch> entries;
  for (const auto &candidate : expected) {
    entries.push_back(StringMatch{
        candidate,
        utils::editDistance(given, candidate)});
  }
  std::sort(entries.begin(), entries.end());
  return entries;
}

bool isKebabStyle(std::string_view sv);

} // namespace precice::utils
