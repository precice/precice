#pragma once

#include <regex>

#include "utils/Event.hpp"

#include <boost/filesystem.hpp>


namespace precice {
namespace utils {

void terminationSignalHandler(int signal)
{
  // Print the events statistics
  precice::utils::EventRegistry::instance().signal_handler(signal);

  // Delete stale .A-B.adress files
  std::regex addrFile(R"(\..*-.*\.address)");

  namespace fs = boost::filesystem;
  fs::path cwd = fs::current_path();

  for (fs::directory_iterator dir(cwd); dir != fs::directory_iterator(); ++dir) {
    std::string file = fs::path(dir->path()).filename().string();
    if (std::regex_match(file, addrFile)) {
      fs::remove(dir->path());
    }
  }

  std::exit(-1);
}
}
}
