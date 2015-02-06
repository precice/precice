// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_ACTION_SHAREDPOINTER_HPP_
#define PRECICE_ACTION_SHAREDPOINTER_HPP_

#include "boost/smart_ptr.hpp"

namespace precice {
namespace action {

class Action;
typedef boost::shared_ptr<Action> PtrAction;

class ActionConfiguration;
typedef boost::shared_ptr<ActionConfiguration> PtrActionConfiguration;

}} // namespace precice, action


#endif /* PRECICE_ACTION_SHAREDPOINTER_HPP_ */
