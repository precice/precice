// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_CONFIG_SHAREDPOINTER_HPP_
#define PRECICE_CONFIG_SHAREDPOINTER_HPP_

#include "boost/smart_ptr.hpp"

namespace precice {
namespace config {

class ParticipantConfiguration;
typedef boost::shared_ptr<ParticipantConfiguration> PtrParticipantConfiguration;

}} // namespace precice, config

#endif /* PRECICE_CONFIG_SHAREDPOINTER_HPP_ */
