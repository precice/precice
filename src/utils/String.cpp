// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "String.hpp"
#include "tarch/Assertions.h"
#include "boost/algorithm/string/split.hpp"
#include "boost/algorithm/string/classification.hpp"

namespace precice {
namespace utils {

std::string wrapText (
  const std::string& text,
  int                linewidth,
  int                indentation )
{
  std::vector<std::string> tokens;
  boost::algorithm::split(tokens, text, boost::algorithm::is_space());
  assertion((int)tokens.size() > 0);
  std::string wrapped;
  int length = 0;
  while (text[length] == ' '){
    wrapped += " ";
    length++;
  }
  for (int i=0; i < (int)tokens.size()-1; i++){
    length += (int)tokens[i].length();
    wrapped += tokens[i];
    if (length + (int)tokens[i+1].length() + 1 > linewidth){
      wrapped += "\n";
      for (int ws=0; ws < indentation; ws++){
        wrapped += " ";
      }
      length = indentation;
    }
    else {
      wrapped += " ";
      length += 1;
    }
  }
  wrapped += tokens.back();
  return wrapped;
}

std::string& checkAppendExtension
(
  std::string&       filename,
  const std::string& extension )
{
  size_t pos = filename.rfind(extension);
  if ((pos == std::string::npos) || (pos != filename.size()-extension.size())){
    filename += extension;
  }
  return filename;
}

}} // namespace precice, utils
