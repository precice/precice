// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at
// http://www5.in.tum.de/wiki/index.php/PreCICE_License

#ifndef PRECICE_NO_MPI

#include "MPIPortsCommunicationFactory.hpp"

#include "MPIPortsCommunication.hpp"

namespace precice {
namespace com {
MPIPortsCommunicationFactory::MPIPortsCommunicationFactory(
    std::string const& addressDirectory)
    : _addressDirectory(addressDirectory) {
  if (_addressDirectory.empty()) {
    _addressDirectory = ".";
  }
}

Communication::SharedPointer
MPIPortsCommunicationFactory::newCommunication() {
  return Communication::SharedPointer(
      new MPIPortsCommunication(_addressDirectory));
}

std::string
MPIPortsCommunicationFactory::addressDirectory() {
  return _addressDirectory;
}
}
} // namespace precice, com

#endif // not PRECICE_NO_MPI
