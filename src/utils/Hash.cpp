#include <boost/uuid/name_generator.hpp>
#include <boost/uuid/string_generator.hpp>
#include <boost/uuid/uuid_io.hpp>

#include "utils/Hash.hpp"
#include "utils/assertion.hpp"

namespace precice::utils {

std::string preciceHash(std::string_view s)
try {
  boost::uuids::string_generator ns_gen;
  auto                           ns = ns_gen("af7ce8f2-a9ee-46cb-38ee-71c318aa3580"); // md5 hash of precice.org as namespace

  boost::uuids::name_generator gen{ns};
  std::string                  hash = boost::uuids::to_string(gen(s.data(), s.size()));
  hash.erase(std::remove(hash.begin(), hash.end(), '-'), hash.end());
  return hash;

} catch (const std::runtime_error &e) {
  PRECICE_UNREACHABLE("preCICE hashing failed", e.what());
  return "";
}

} // namespace precice::utils
