// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at
// http://www5.in.tum.de/wiki/index.php/PreCICE_License

#ifndef PRECICE_NO_MPI

#ifndef PRECICE_COM_MPI_PORTS_COMMUNICATION_FACTORY_HPP_
#define PRECICE_COM_MPI_PORTS_COMMUNICATION_FACTORY_HPP_

#include "CommunicationFactory.hpp"

#include <string>

namespace precice {
namespace com {
class MPIPortsCommunicationFactory : public CommunicationFactory {
public:
  /**
   * @brief Constructor.
   */
  MPIPortsCommunicationFactory(std::string const& addressDirectory = ".");

  Communication::SharedPointer newCommunication();

  std::string addressDirectory();

private:
  std::string _addressDirectory;
};
}
} // namespace precice, com

#endif /* PRECICE_COM_MPI_PORTS_COMMUNICATION_FACTORY_HPP_ */

#endif // not PRECICE_NO_MPI
