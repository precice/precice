// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_COM_SHAREDPOINTER_HPP_
#define PRECICE_COM_SHAREDPOINTER_HPP_

#include "boost/smart_ptr.hpp"

namespace precice {
namespace com {

class Communication;
typedef boost::shared_ptr<Communication> PtrCommunication;

class CommunicationConfiguration;
typedef boost::shared_ptr<CommunicationConfiguration> PtrCommunicationConfiguration;

}} // namespace precice, com

#endif /* SHAREDPOINTER_HPP_ */
