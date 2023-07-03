#include <algorithm>
#include <boost/algorithm/string/trim.hpp>
#include <boost/filesystem.hpp>
#include <boost/uuid/name_generator.hpp>
#include <boost/uuid/string_generator.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <chrono>
#include <fstream>
#include <stdexcept>
#include <thread>

#include "com/ConnectionInfoPublisher.hpp"
#include "logging/LogMacros.hpp"
#include "precice/types.hpp"
#include "utils/assertion.hpp"

namespace bfs = boost::filesystem;

namespace precice::com {

namespace {

std::string preciceFancyHash(const std::string &s)
try {
  boost::uuids::string_generator ns_gen;
  auto                           ns = ns_gen("af7ce8f2-a9ee-46cb-38ee-71c318aa3580"); // md5 hash of precice.org as namespace

  boost::uuids::name_generator gen{ns};
  return boost::uuids::to_string(gen(s));

} catch (const std::runtime_error &e) {
  PRECICE_UNREACHABLE("preCICE hashing failed", e.what());
  return "";
}
} // namespace

std::string impl::hashedFilePath(const std::string &acceptorName, const std::string &requesterName, const std::string &tag, Rank rank)
{
  constexpr int     firstLevelLen = 2;
  std::string const s             = acceptorName + tag + requesterName + std::to_string(rank);
  std::string       hash          = preciceFancyHash(s);
  hash.erase(std::remove(hash.begin(), hash.end(), '-'), hash.end());

  auto p = bfs::path(hash.substr(0, firstLevelLen)) / hash.substr(firstLevelLen);

  return p.string();
}

std::string impl::localDirectory(const std::string &acceptorName, const std::string &requesterName, const std::string &addressDirectory)
{
  std::string directional = acceptorName + "-" + requesterName;

  auto p = bfs::path(addressDirectory) / "precice-run" / directional;

  return p.string();
}

std::string ConnectionInfoPublisher::getLocalDirectory() const
{
  return impl::localDirectory(acceptorName, requesterName, addressDirectory);
}

std::string ConnectionInfoPublisher::getFilename() const
{
  auto local  = getLocalDirectory();
  auto hashed = impl::hashedFilePath(acceptorName, requesterName, tag, rank);
  auto p      = bfs::path(getLocalDirectory()) / hashed;

  return p.string();
}

std::string ConnectionInfoReader::read() const
{
  auto path = getFilename();

  PRECICE_DEBUG("Waiting for connection file \"{}\"", path);
  const auto waitdelay = std::chrono::milliseconds(1);
  while (!bfs::exists(path)) {
    std::this_thread::sleep_for(waitdelay);
  }
  PRECICE_ASSERT(bfs::exists(path));
  PRECICE_DEBUG("Found connection file \"{}\"", path);

  std::ifstream ifs(path);
  PRECICE_CHECK(ifs,
                "Unable to establish connection as the connection file \"{}\" couldn't be opened.",
                path);
  std::string addressData;
  std::getline(ifs, addressData);
  PRECICE_CHECK(!addressData.empty(),
                "Unable to establish connection as the connection file \"{}\" is empty. "
                "Please report this bug to the preCICE developers.",
                path);
  boost::algorithm::trim_right(addressData);
  return addressData;
}

ConnectionInfoWriter::~ConnectionInfoWriter()
{
  bfs::path path(getFilename());
  if (!bfs::exists(path)) {
    PRECICE_WARN("Cannot clean-up the connection file \"{}\" as it doesn't exist. "
                 "In case of connection problems, please report this to the preCICE developers.",
                 path.generic_string());
    return;
  }
  PRECICE_DEBUG("Deleting connection file \"{}\"", path.generic_string());
  try {
    bfs::remove(path);
    if (bfs::exists(path)) {
      PRECICE_WARN("The connection file \"{}\" wasn't properly removed. "
                   "Make sure to delete the \"precice-run\" directory before restarting the simulation.",
                   path.generic_string());
    }
  } catch (const bfs::filesystem_error &e) {
    PRECICE_WARN("Unable to clean-up connection file due to error: {}. "
                 "Make sure to delete the \"precice-run\" directory before restarting the simulation.",
                 e.what());
  }
}

void ConnectionInfoWriter::write(std::string const &info) const
{
  auto path = getFilename();
  auto tmp  = bfs::path(path + "~");

  {
    auto message = "Unable to establish connection as a {}connection file already exists at \"{}\". "
                   "This is likely a leftover of a previous crash or stop during communication build-up. "
                   "Please remove the \"precice-run\" directory and restart the simulation.";
    PRECICE_CHECK(!bfs::exists(path), message, "", path);
    PRECICE_CHECK(!bfs::exists(tmp), message, "temporary ")
  }

  PRECICE_DEBUG("Writing temporary connection file \"{}\"", tmp.generic_string());
  bfs::create_directories(tmp.parent_path());
  {
    std::ofstream ofs(tmp.string());
    PRECICE_CHECK(ofs, "Unable to establish connection as the temporary connection file \"{}\" couldn't be opened.", tmp.generic_string());
    fmt::print(ofs,
               "{}\nAcceptor: {}, Requester: {}, Tag: {}, Rank: {}",
               info, acceptorName, requesterName, tag, rank);
  }
  PRECICE_CHECK(bfs::exists(tmp),
                "Unable to establish connection as the temporary connection file \"{}\" was written, but doesn't exist on disk. "
                "Please report this bug to the preCICE developers.",
                tmp.generic_string());

  PRECICE_DEBUG("Publishing connection file \"{}\"", path);
  bfs::rename(tmp, path);
  if (bfs::exists(tmp)) {
    PRECICE_WARN("The temporary connection file \"{}\" wasn't properly removed. "
                 "Make sure to delete the \"precice-run\" directory before restarting the simulation.",
                 tmp.generic_string());
  }
  PRECICE_CHECK(bfs::exists(path),
                "Unable to establish connection as the connection file \"{}\" doesn't exist on disk. "
                "Please report this bug to the preCICE developers.",
                path);
}

} // namespace precice::com
