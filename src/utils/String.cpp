// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "String.hpp"
#include "tarch/Assertions.h"

namespace precice {
namespace utils {

std::string wrapText (
  const std::string& text,
  int                linewidth,
  int                indentation )
{
  std::vector<std::string> tokens = tokenize(text, " ");
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

std::vector<std::string> tokenize
(
   const std::string & toTokenize,
   const std::string & delimiter )
{
   std::string::size_type lengthDelimiter = delimiter.size();
   std::vector<std::string> tokens;
   std::string::size_type startPosition = 0;
   std::string::size_type location = toTokenize.find ( delimiter, startPosition );
   while ( location != std::string::npos ) {
      if ( location > startPosition ) {
         tokens.push_back ( toTokenize.substr(startPosition, location - startPosition) );
      }
      startPosition = location + lengthDelimiter;
      location = toTokenize.find ( delimiter, startPosition );
   }
   if ( (location > startPosition) && (startPosition != toTokenize.size()) ) {
      tokens.push_back ( toTokenize.substr(startPosition, location - startPosition) );
   }
   return tokens;
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
