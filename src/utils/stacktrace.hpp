#pragma once

#include <execinfo.h>	// for backtrace
#include <cxxabi.h>	// for __cxa_demangle

#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>

/// Prints a demangled stack backtrace of the caller function
static inline std::string getStacktrace()
{
  std::ostringstream strm;
  const int max_frames = 100;
  strm << "Stack Trace:" << std::endl;

  // Storage array for stack trace address data
  void* addrlist[max_frames+1];

  // Retrieve current stack addresses
  int addrlen = backtrace(addrlist, sizeof(addrlist) / sizeof(void*));

  if (addrlen == 0) {
    strm << "  <empty, possibly corrupt>" << std::endl;
    return strm.str();
  }

  // Resolve addresses into strings containing "filename(function+address)", this array must be free()-ed
  char** symbollist = backtrace_symbols(addrlist, addrlen);

  // Allocate string which will be filled with the demangled function name
  size_t funcnamesize = 256;
  char* funcname = new char[funcnamesize];
  
  // Iterate over the returned symbol lines. skip the first, it is the address of this function.
  for (int i = 1; i < addrlen; i++) {
    std::string strSymbols(symbollist[i]);
    auto op = strSymbols.find("(");
    auto cp = strSymbols.find(")");
    auto plus = strSymbols.find("+");
    std::string strSymbol = strSymbols.substr(op+1, plus-op-1);
    std::string strOffset = strSymbols.substr(plus+1, cp-plus-1);

    int status;
    char* ret = abi::__cxa_demangle(strSymbol.c_str(), funcname, &funcnamesize, &status);
    strm << "  (" << std::setw(2) << addrlen-i << ")";
    if (status == 0) {
      funcname = ret; // Use possibly realloc()-ed string
      strm << "  " << strSymbols.substr(0, op) << " : " << funcname << "+" << strOffset <<  std::endl;
    }
    else {
      // Demangling failed. Output function name as a C function with no arguments.
      strm << "  " << strSymbols.substr(0, op) << " : " << strSymbol << "()+" << strOffset <<  std::endl;
    }
  }

  delete[] funcname;
  free(symbollist);
  
  return strm.str();
}
