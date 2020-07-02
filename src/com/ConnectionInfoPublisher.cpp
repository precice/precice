#include "ConnectionInfoPublisher.hpp"
#include <algorithm>
#include <boost/filesystem.hpp>
#include <boost/uuid/name_generator.hpp>
#include <boost/uuid/string_generator.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <chrono>
#include <istream>
#include <thread>
#include "logging/LogMacros.hpp"

namespace precice {
namespace com {

std::string impl::hashedFilePath(const std::string &acceptorName, const std::string &requesterName, const std::string &tag, int rank)
{
  using namespace boost::filesystem;

  constexpr int                  firstLevelLen = 2;
  boost::uuids::string_generator ns_gen;
  auto                           ns = ns_gen("af7ce8f2-a9ee-46cb-38ee-71c318aa3580"); // md5 hash of precice.org as namespace

  boost::uuids::name_generator gen{ns};
  std::string const            s    = acceptorName + tag + requesterName + std::to_string(rank);
  std::string                  hash = boost::uuids::to_string(gen(s));
  hash.erase(std::remove(hash.begin(), hash.end(), '-'), hash.end());

  path p = path(hash.substr(0, firstLevelLen)) / hash.substr(firstLevelLen);

  return p.string();
}

std::string impl::localDirectory(const std::string &acceptorName, const std::string &requesterName, const std::string &addressDirectory)
{
  using namespace boost::filesystem;
  std::string directional = acceptorName + "-" + requesterName;

  path p = path(addressDirectory) / path("precice-run") / path(directional);

  return p.string();
}

std::string ConnectionInfoPublisher::getLocalDirectory() const
{
  return impl::localDirectory(acceptorName, requesterName, addressDirectory);
}

std::string ConnectionInfoPublisher::getFilename() const
{
  using namespace boost::filesystem;

  auto local  = getLocalDirectory();
  auto hashed = impl::hashedFilePath(acceptorName, requesterName, tag, rank);
  path p      = path(getLocalDirectory()) / path(hashed);

  return p.string();
}

std::string ConnectionInfoReader::read() const
{
  std::ifstream ifs;
  auto          path = getFilename();
  PRECICE_DEBUG("Waiting for connection file " << path);
  do {
    ifs.open(path, std::ifstream::in);
    std::this_thread::sleep_for(std::chrono::milliseconds(1));
  } while (not ifs);
  PRECICE_DEBUG("Found connection file " << path);

  std::string addressData;
  ifs >> addressData;
  return addressData;
}

ConnectionInfoWriter::~ConnectionInfoWriter()
{
  namespace fs = boost::filesystem;
  fs::path p(getFilename());
  PRECICE_DEBUG("Deleting connection file " << p.string());
  fs::remove(p);
}

void ConnectionInfoWriter::write(std::string const &info) const
{
  namespace fs = boost::filesystem;
  auto path    = getFilename();
  auto tmp     = fs::path(path + "~");
  PRECICE_DEBUG("Writing connection file " << path);
  fs::create_directories(tmp.parent_path());
  {
    std::ofstream ofs(tmp.string(), std::ofstream::out);
    ofs << info << "\n";
    ofs << "Acceptor: " << acceptorName << ", Requester: " << requesterName << ", Tag: " << tag << ", Rank: " << rank << "\n";
  }
  fs::rename(tmp, path);
}

} // namespace com
} // namespace precice
