#pragma once

#include <string>

namespace precice {
namespace utils {

std::string wrapText (
  const std::string& text,
  int                linewidth,
  int                indentation );

/**
 * @brief Checks if filename has the given extension, if not appends it.
 *
 * @return filename with extension.
 */
std::string & checkAppendExtension(std::string& filename, const std::string& extension);

}} // namespace precice, utils

