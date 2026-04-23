#include "CommunicationConfiguration.hpp"
#include <memory>
#include <ostream>
#include "com/MPIDirectCommunication.hpp"
#include "com/MPIIntraComm.hpp"
#include "com/MPIPortsCommunication.hpp"
#include "com/SocketCommunication.hpp"
#include "com/SocketIntraComm.hpp"
#include "logging/LogMacros.hpp"
#include "utils/Helpers.hpp"
#include "utils/assertion.hpp"
#include "xml/XMLTag.hpp"

namespace precice::com {
PtrIntraCommunication CommunicationConfiguration::createIntraCommunication(
    const xml::XMLTag &tag) const
{
  PtrIntraCommunication comm;
  if (tag.getName() == "sockets") {
    std::string network = tag.getStringAttributeValue("network");
    int         port    = tag.getIntAttributeValue("port");

    PRECICE_CHECK(utils::isValidPort(port),
                  "A sockets communication was configured with an invalid port \"{}\". "
                  "Please check the \"ports=\" attributes of your socket connections.",
                  port);

    std::string dir = tag.getStringAttributeValue("exchange-directory");
    comm            = std::make_shared<com::SocketIntraComm>(port, false, network, dir);
  } else if (tag.getName() == "mpi") {
#ifdef PRECICE_NO_MPI
    PRECICE_ERROR("Communication type \"mpi\" can only be used if preCICE was compiled with MPI support enabled. "
                  "Either switch to a \"sockets\" communication or recompile preCICE with \"PRECICE_MPICommunication=ON\".");
#else
    comm = std::make_shared<com::MPIIntraComm>();
#endif
  }
  PRECICE_ASSERT(comm != nullptr);
  return comm;
}

} // namespace precice::com
