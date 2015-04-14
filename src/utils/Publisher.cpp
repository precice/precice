#include "Publisher.hpp"

#include "EventTimings.hpp"

#include <fstream>
#include <sstream>

namespace precice {
namespace utils {
std::string Publisher::_prefix;

Publisher::ScopedPrefix::ScopedPrefix(std::string const& prefix)
    : _prefix(Publisher::prefix()) {
  Publisher::setPrefix(prefix);
}

Publisher::ScopedPrefix::~ScopedPrefix() {
  Publisher::setPrefix(_prefix);
}

Publisher::ScopedPublication::ScopedPublication(std::string const& filePath,
                                                std::string const& data)
    : _filePath(filePath) {
  Publisher::write(_filePath, data);
}

Publisher::ScopedPublication::~ScopedPublication() {
  Publisher::remove(_filePath);
}

void
Publisher::read(std::string const& filePath, std::string& data) {
  std::ifstream ifs;

  auto beforeReadTimeStamp =
      std::chrono::duration_cast<std::chrono::milliseconds>(
          Event::Clock::now().time_since_epoch());

  do {
    ifs.open(filePath, std::ifstream::in);
  } while (not ifs);

  auto afterReadTimeStamp =
      std::chrono::duration_cast<std::chrono::milliseconds>(
          Event::Clock::now().time_since_epoch());

  std::chrono::milliseconds::rep writeTimeStampCount;

  {
    std::string line;

    std::getline(ifs, line);

    std::istringstream iss(line);

    iss >> writeTimeStampCount;

    std::getline(ifs, data);
  }

  std::chrono::milliseconds writeTimeStamp(writeTimeStampCount);

  std::chrono::milliseconds readDuration;

  if (writeTimeStamp > beforeReadTimeStamp)
    readDuration = afterReadTimeStamp - writeTimeStamp;
  else
    readDuration = afterReadTimeStamp - beforeReadTimeStamp;

  std::string eventName = _prefix;

  if (not _prefix.empty())
    eventName += "/";

  eventName += "Publisher::read";

  Event(eventName, readDuration);
}

void
Publisher::write(std::string const& filePath, std::string const& data) {
  {
    std::ofstream ofs(filePath + "~", std::ofstream::out);

    auto writeTimeStamp = std::chrono::duration_cast<std::chrono::milliseconds>(
        Event::Clock::now().time_since_epoch());

    ofs << writeTimeStamp.count() << "\n" << data;
  }

  std::rename((filePath + "~").c_str(), filePath.c_str());
}

void
Publisher::remove(std::string const& filePath) {
  std::remove(filePath.c_str());
}

void
Publisher::setPrefix(std::string const& prefix) {
  _prefix = prefix;
}

std::string const&
Publisher::prefix() {
  return _prefix;
}
} // namespace utils
} // namespace precice
