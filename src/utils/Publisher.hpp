// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at
// http://www5.in.tum.de/wiki/index.php/PreCICE_License

#ifndef PRECICE_UTILS_PUBLISHER_HPP_
#define PRECICE_UTILS_PUBLISHER_HPP_

#include "Helpers.hpp"

#include <string>

namespace precice {
namespace utils {
class Publisher {
public:
  struct ScopedSetEventNamePrefix {
    ScopedSetEventNamePrefix(std::string const& prefix);

    ~ScopedSetEventNamePrefix();

  private:
    std::string _prefix;
  };

  struct ScopedPushDirectory {
    ScopedPushDirectory(std::string const& dp);

    ~ScopedPushDirectory();
  };

  struct ScopedChangePrefixDirectory {
    ScopedChangePrefixDirectory(std::string const& pdp);

    ~ScopedChangePrefixDirectory();

  private:
    std::string _pdp;
  };

public:
  static std::string parentPath(std::string const& p);

  static bool createDirectory(std::string const& dp);

  static bool exists(std::string const& p);

  static bool remove(std::string const& p);

  static void rename(std::string const& op, std::string const& np);

  static bool pushDirectory(std::string const& dp);

  static bool popDirectory();

  static void changePrefixDirectory(std::string const& pdp);

  static std::string const& prefixDirectoryPath();

  static void setEventNamePrefix(std::string const& prefix);

  static std::string const& eventNamePrefix();

public:
  Publisher(std::string const& fp);

  void read(std::string& data) const;

  void write(std::string const& data) const;

  std::string const& filePath() const;

private:
  static std::string buildFilePath(std::string const& fp);

private:
  static std::string _pdp;

  static Stack<std::string> _dps;

  static std::string _prefix;

private:
  std::string _fp;
};

class ScopedPublisher : public Publisher {
public:
  ScopedPublisher(std::string const& fp);

  ~ScopedPublisher();
};

} // namespace utils
} // namespace precice

#endif /* PRECICE_UTILS_PUBLISHER_HPP_ */
