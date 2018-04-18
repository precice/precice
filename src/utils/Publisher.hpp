#pragma once

#include <stack>
#include <string>

namespace precice
{
namespace utils
{

/**
 * @brief Publisher Class. This utility class can be used to publish connection information and
 * was established to separate technical details from the communication.
 * Sub-folders can be used to allow for a more efficient reading and writing.
 *
 * @todo implementation should be substituted by a proper publishing strategy, not via files.
 * (to allow for exa-scale one day)
 */
class Publisher
{
public:
  struct ScopedSetEventNamePrefix {
    explicit ScopedSetEventNamePrefix(std::string const &prefix);

    ~ScopedSetEventNamePrefix();

  private:
    std::string _prefix;
  };

  struct ScopedPushDirectory {
    explicit ScopedPushDirectory(std::string const &dp);

    ~ScopedPushDirectory();
  };

  /// Sets the address publishing directory. Resets it, when object leaves scope.
  struct ScopedChangePrefixDirectory {
    explicit ScopedChangePrefixDirectory(std::string const &pdp);

    ~ScopedChangePrefixDirectory();

    std::string _pdp;
  };

public:
  static std::string parentPath(std::string const &p);

  static bool createDirectory(std::string const &dp);

  static bool exists(std::string const &p);

  static bool remove(std::string const &p);

  static void rename(std::string const &op, std::string const &np);

  static bool pushDirectory(std::string const &dp);

  static bool popDirectory();

  static void changePrefixDirectory(std::string const &pdp);

  static std::string const &prefixDirectoryPath();

  static void setEventNamePrefix(std::string const &prefix);

  static std::string const &eventNamePrefix();

public:
  explicit Publisher(std::string const &fp);

  std::string read() const;

  void write(std::string const &data) const;

  std::string const &filePath() const;

private:
  static std::string buildFilePath(std::string const &fp);

  static std::string _pdp;

  static std::stack<std::string> _dps;

  static std::string _prefix;

private:
  std::string _fp;
};

class ScopedPublisher : public Publisher
{
public:
  explicit ScopedPublisher(std::string const &fp);

  ~ScopedPublisher();
};

} // namespace utils
} // namespace precice
