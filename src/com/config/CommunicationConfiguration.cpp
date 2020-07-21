#include "CommunicationConfiguration.hpp"
#include <memory>
#include <ostream>
#include "com/MPIDirectCommunication.hpp"
#include "com/MPIPortsCommunication.hpp"
#include "com/SocketCommunication.hpp"
#include "logging/LogMacros.hpp"
#include "utils/Helpers.hpp"
#include "utils/assertion.hpp"
#include "xml/XMLTag.hpp"

namespace precice {
namespace com {
PtrCommunication CommunicationConfiguration::createCommunication(
    const xml::XMLTag &tag) const
{
  com::PtrCommunication com;
  if (tag.getName() == "sockets") {
    std::string network = tag.getStringAttributeValue("network");
    int         port    = tag.getIntAttributeValue("port");

    PRECICE_CHECK(utils::isValidPort(port),
                  "A sockets communication was configured with an invalid port \"" << port << "\". Please check the \"ports="
                                                                                              "\" attributes of your socket connections.");

    std::string dir = tag.getStringAttributeValue("exchange-directory");
    com             = std::make_shared<com::SocketCommunication>(port, false, network, dir);
  } else if (tag.getName() == "mpi") {
    std::string dir = tag.getStringAttributeValue("exchange-directory");
#ifdef PRECICE_NO_MPI
    PRECICE_ERROR("Communication type \"mpi\" can only be used if preCICE was compiled with MPI support enabled. "
                  "Either switch to a \"sockets\" communication or recompile preCICE with \"PRECICE_MPICommunication=ON\".");
#else
    com = std::make_shared<com::MPIPortsCommunication>(dir);
#endif
  } else if (tag.getName() == "mpi-single") {
#ifdef PRECICE_NO_MPI
    PRECICE_ERROR("Communication type \"mpi-single\" can only be used if preCICE was compiled with MPI support enabled. "
                  "Either switch to a \"sockets\" communication or recompile preCICE with \"PRECICE_MPICommunication=ON\".");
#else
    com = std::make_shared<com::MPIDirectCommunication>();
#endif
  }
  PRECICE_ASSERT(com != nullptr);
  return com;
}

} // namespace com
} // namespace precice
