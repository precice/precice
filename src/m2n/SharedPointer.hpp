// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License

#ifndef PRECICE_M2N_SHAREDPOINTER_HPP_
#define PRECICE_M2N_SHAREDPOINTER_HPP_

#include <memory>

namespace precice {
namespace m2n {
class DistributedCommunication;
class DistributedComFactory;
class M2N;
class M2NConfiguration;

using PtrDistributedCommunication = std::shared_ptr<DistributedCommunication>;
using PtrDistributedComFactory    = std::shared_ptr<DistributedComFactory>;
using PtrM2N                      = std::shared_ptr<M2N>;
using PtrM2NConfiguration         = std::shared_ptr<M2NConfiguration>;
}} // namespace precice, m2n

#endif /* PRECICE_M2N_SHAREDPOINTER_HPP_ */
