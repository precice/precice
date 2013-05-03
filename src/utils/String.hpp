// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_UTILS_STRING_HPP_
#define PRECICE_UTILS_STRING_HPP_

#include <vector>
#include <string>

namespace precice {
namespace utils {

std::string wrapText (
  const std::string& text,
  int                linewidth,
  int                indentation );

/**
 * @brief Tokenizes a string into substrings separated by delimiters.
 *
 * An example:
 * - string to tokenize    : "example++string++to++tokenize"
 * - delimiter             : "++"
 * - result vector entries : "example", "string", "to", "tokenize"
 *
 * @param toTokenize [IN] String to be tokenized.
 * @param delimiter  [IN] Delimiting characters, used to separate tokens.
 */
std::vector<std::string> tokenize (
  const std::string & toTokenize,
  const std::string & delimiter );

/**
 * @brief Checks if filename has the given extension, if not appends it.
 *
 * @return filename with extension.
 */
std::string& checkAppendExtension (
  std::string& filename,
  const std::string& extension );

}} // namespace precice, utils

#endif /* PRECICE_UTILS_STRING_HPP_ */
