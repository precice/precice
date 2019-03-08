#include "ConnectionInfoPublisher.hpp"
#include <chrono>
#include <thread>
#include <boost/uuid/name_generator.hpp>
#include <boost/uuid/string_generator.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <boost/filesystem.hpp>

namespace precice
{
namespace com
{

std::string ConnectionInfoPublisher::getFilename() const
{
  using namespace boost::filesystem;
  
  constexpr int firstLevelLen = 2;
  boost::uuids::string_generator ns_gen;
  auto ns = ns_gen("af7ce8f2-a9ee-46cb-38ee-71c318aa3580"); // md5 hash of precice.org as namespace
  
  boost::uuids::name_generator gen{ns};
  std::string const s = acceptorName + requesterName + std::to_string(rank);
  std::string hash = boost::uuids::to_string(gen(s));
  hash.erase(std::remove(hash.begin(), hash.end(), '-'), hash.end());
  
  path p = path(addressDirectory)
    / path(".precice")
    / path(hash.substr(0, firstLevelLen))
    / hash.substr(firstLevelLen);
  
  return p.string();
}


std::string ConnectionInfoReader::read() const
{
  std::ifstream ifs;
  auto path = getFilename();
  do {
    ifs.open(path, std::ifstream::in);
    std::this_thread::sleep_for(std::chrono::milliseconds(1));
  } while (not ifs);
  
  std::string addressData;
  ifs >> addressData;
  return addressData;
}


ConnectionInfoWriter::~ConnectionInfoWriter()
{
  namespace fs = boost::filesystem;
  fs::path p(getFilename());
  try {
    fs::remove(getFilename());
    if (fs::is_empty(p.parent_path()))
      fs::remove(p.parent_path()); // also remove parent dir, i.e, first part of hash
    if (fs::is_empty(p.parent_path().parent_path()))
      fs::remove(p.parent_path().parent_path()); // and also the .precice, only if empty
  }
  catch (fs::filesystem_error const & e) {
    WARN("Filesystem error when deleting connection info files: " << e.what());
  }
}


void ConnectionInfoWriter::write(std::string const & info) const
{
  namespace fs = boost::filesystem;
  auto path = getFilename();
  fs::create_directories(fs::path(path).parent_path());
  std::ofstream ofs(path + "~");
  ofs << info;
  ofs.close();
  fs::rename(path + "~", path);
}

}
}
