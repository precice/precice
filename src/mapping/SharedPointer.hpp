// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_MAPPING_SHAREDPOINTER_HPP_
#define PRECICE_MAPPING_SHAREDPOINTER_HPP_

#include "boost/smart_ptr.hpp"

namespace precice {
namespace mapping {

class Mapping;
typedef boost::shared_ptr<Mapping> PtrMapping;

class MappingConfiguration;
typedef boost::shared_ptr<MappingConfiguration> PtrMappingConfiguration;

}} // namsepace precice, mapping

#endif /* PRECICE_MAPPING_SHAREDPOINTER_HPP_ */
