// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_M2N_SHAREDPOINTER_HPP_
#define PRECICE_M2N_SHAREDPOINTER_HPP_

#include "boost/smart_ptr.hpp"

namespace precice {
namespace m2n {

class DistributedCommunication;
typedef boost::shared_ptr<DistributedCommunication> PtrDistributedCommunication;

class DistributedComFactory;
typedef boost::shared_ptr<DistributedComFactory> PtrDistributedComFactory;

class M2N;
typedef boost::shared_ptr<M2N> PtrM2N;

class M2NConfiguration;
typedef boost::shared_ptr<M2NConfiguration> PtrM2NConfiguration;

}} // namespace precice, m2n

#endif /* PRECICE_M2N_SHAREDPOINTER_HPP_ */
