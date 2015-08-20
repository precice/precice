#pragma once

#include <iostream>
using namespace std;
#include <regex>

#include "utils/EventTimings.hpp"
#include "utils/Globals.hpp"

#include <boost/filesystem.hpp>


namespace precice {
namespace utils {

tarch::logging::Log _log("precice::utils::terminationSignalHandler");


void terminationSignalHandler(int signal)
{
  cout << "SIGNAL HANDLER CALLED" << endl;
  precice::utils::EventRegistry::signal_handler(signal);

  boost::filesystem::path cwd = boost::filesystem::current_path();

  std::regex addrFile(R"(\..*-.*\.address)");
  for (boost::filesystem::directory_entry& dirEntry : boost::filesystem::directory_iterator(cwd)) {
    std::string file = boost::filesystem::path(dirEntry).filename().string();
    if (std::regex_match(file, addrFile)) {
      boost::filesystem::remove(dirEntry);
    }
  }
}
}
}
