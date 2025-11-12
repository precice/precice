#include <algorithm>
#include <boost/algorithm/string/trim.hpp>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <thread>

#include "com/ConnectionInfoPublisher.hpp"
#include "logging/LogMacros.hpp"
#include "precice/impl/Types.hpp"
#include "utils/Hash.hpp"
#include "utils/assertion.hpp"

#include <profiling/Event.hpp>

using precice::profiling::Event;

namespace fs = std::filesystem;
namespace precice::com {

std::string impl::hashedFilePath(std::string_view acceptorName, std::string_view requesterName, std::string_view tag, Rank rank)
{
  constexpr int     firstLevelLen = 2;
  std::string const s             = std::string(acceptorName).append(tag).append(requesterName).append(std::to_string(rank));
  std::string       hash          = utils::preciceHash(s);

  auto p = fs::path(hash.substr(0, firstLevelLen)) / hash.substr(firstLevelLen);

  return p.string();
}

std::string impl::localDirectory(std::string_view acceptorName, std::string_view requesterName, std::string_view addressDirectory)
{
  std::string directional = std::string(acceptorName).append("-").append(requesterName);

  auto p = fs::path(addressDirectory.begin(), addressDirectory.end()) / "precice-run" / directional;

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
  auto p      = fs::path(getLocalDirectory()) / hashed;

  return p.string();
}

std::string ConnectionInfoReader::read() const
{
  auto path = getFilename();

  Event e0 = Event("ConnectionInfoReader.read.waitingForConnectionFile");
  PRECICE_DEBUG("Waiting for connection file \"{}\"", path);
  const auto waitdelay = std::chrono::milliseconds(1);
  while (!fs::exists(path)) {
    std::this_thread::sleep_for(waitdelay);
  }
  PRECICE_ASSERT(fs::exists(path));
  PRECICE_DEBUG("Found connection file \"{}\"", path);

  e0.stop();
  Event e1 = Event("ConnectionInfoReader.read.readingConnectionFile");

  std::ifstream ifs(path);

  if (!ifs) {
    PRECICE_DEBUG("Opening connection file \"{}\" failed. Retrying", path);
    std::this_thread::sleep_for(waitdelay);
    ifs.clear();
    ifs.open(path);
  }

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

  e1.stop();

  return addressData;
}

ConnectionInfoWriter::~ConnectionInfoWriter()
{
  Event e0 = Event("ConnectionInfoWriter.init.constructingFilename");

  fs::path path(getFilename());

  e0.stop();
  Event e1 = Event("ConnectionInfoWriter.init.checkingConnectingFile");

  if (!fs::exists(path)) {
    PRECICE_WARN("Cannot clean-up the connection file \"{}\" as it doesn't exist. "
                 "In case of connection problems, please report this to the preCICE developers.",
                 path.generic_string());
    return;
  }

  e1.stop();
  Event e2 = Event("ConnectionInfoWriter.init.deletingConnectionFile");

  PRECICE_DEBUG("Deleting connection file \"{}\"", path.generic_string());
  try {
    fs::remove(path);
    PRECICE_WARN_IF(
        fs::exists(path),
        "The connection file \"{}\" wasn't properly removed. "
        "Make sure to delete the \"precice-run\" directory before restarting the simulation.",
        path.generic_string());
  } catch (const fs::filesystem_error &e) {
    PRECICE_WARN("Unable to clean-up connection file due to error: {}. "
                 "Make sure to delete the \"precice-run\" directory before restarting the simulation.",
                 e.what());
  }

  e2.stop();
}

void ConnectionInfoWriter::write(std::string_view info) const
{
  Event e0 = Event("ConnectionInfoWriter.write.constructingFilename");
  auto path = getFilename();
  auto tmp  = fs::path(path + "~");

  e0.stop();
  Event e1 = Event("ConnectionInfoWriter.write.checkExistingFile");

  {
    auto message = "Unable to establish connection as a {}connection file already exists at \"{}\". "
                   "This is likely a leftover of a previous crash or stop during communication build-up. "
                   "Please remove the \"precice-run\" directory and restart the simulation.";
    PRECICE_CHECK(!fs::exists(path), message, "", path);
    PRECICE_CHECK(!fs::exists(tmp), message, "temporary ");
  }

  e1.stop();
  Event e2 = Event("ConnectionInfoWriter.write.writingTemporaryFile");

  PRECICE_DEBUG("Writing temporary connection file \"{}\"", tmp.generic_string());
  fs::create_directories(tmp.parent_path());
  {
    std::ofstream ofs(tmp);

    if (!ofs) {
      PRECICE_DEBUG("Opening temporary connection file \"{}\" failed. Retrying", tmp.generic_string());
      std::this_thread::sleep_for(std::chrono::milliseconds{1});
      ofs.clear();
      ofs.open(tmp);
    }

    PRECICE_CHECK(ofs, "Unable to establish connection as the temporary connection file \"{}\" couldn't be opened.", tmp.generic_string());
    fmt::print(ofs,
               "{}\nAcceptor: {}, Requester: {}, Tag: {}, Rank: {}",
               info, acceptorName, requesterName, tag, rank);
  }
  e2.stop();
  Event e3 = Event("ConnectionInfoWriter.write.checkingTemporaryFile");
  PRECICE_CHECK(fs::exists(tmp),
                "Unable to establish connection as the temporary connection file \"{}\" was written, but doesn't exist on disk. "
                "Please report this bug to the preCICE developers.",
                tmp.generic_string());
  e3.stop();

  Event e4 = Event("ConnectionInfoWriter.write.publishingFile");

  PRECICE_DEBUG("Publishing connection file \"{}\"", path);
  fs::rename(tmp, path);

  e4.stop();
  Event e5 = Event("ConnectionInfoWriter.write.checkingFile");

  PRECICE_WARN_IF(
      fs::exists(tmp),
      "The temporary connection file \"{}\" wasn't properly removed. "
      "Make sure to delete the \"precice-run\" directory before restarting the simulation.",
      tmp.generic_string());
  PRECICE_CHECK(fs::exists(path),
                "Unable to establish connection as the connection file \"{}\" doesn't exist on disk. "
                "Please report this bug to the preCICE developers.",
                path);

  e5.stop();
}

} // namespace precice::com
