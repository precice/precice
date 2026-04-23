#pragma once

#include <memory>

namespace precice::com {

class Communication;
class CommunicationFactory;
class IntraCommunication;
class Request;

using PtrCommunication        = std::shared_ptr<Communication>;
using PtrCommunicationFactory = std::shared_ptr<CommunicationFactory>;
using PtrIntraCommunication   = std::shared_ptr<IntraCommunication>;
using PtrRequest              = std::shared_ptr<Request>;
} // namespace precice::com
