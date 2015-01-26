// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_NEWSPACETREE_SHAREDPOINTER_HPP_
#define PRECICE_NEWSPACETREE_SHAREDPOINTER_HPP_

#include "boost/smart_ptr.hpp"

namespace precice {
namespace spacetree {

class Spacetree;
typedef boost::shared_ptr<Spacetree> PtrSpacetree;

class SpacetreeConfiguration;
typedef boost::shared_ptr<SpacetreeConfiguration> PtrSpacetreeConfiguration;


}} // namespace precice, spacetree


#endif /* PRECICE_NEWSPACETREE_SHAREDPOINTER_HPP_ */
