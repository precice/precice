// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License

#ifndef PRECICE_COM_SHAREDPOINTER_HPP_
#define PRECICE_COM_SHAREDPOINTER_HPP_

#include <memory>

namespace precice {
namespace com {
class Communication;
class CommunicationFactory;
class Request;

using PtrCommunication        = std::shared_ptr<Communication>;
using PtrCommunicationFactory = std::shared_ptr<CommunicationFactory>;
using PtrRequest              = std::shared_ptr<Request>;
}} // namespace precice, com

#endif /* PRECICE_COM_SHAREDPOINTER_HPP_ */
