#pragma once

#include <memory>

namespace precice {
namespace com {

class Communication;
class CommunicationFactory;
class Request;

using PtrCommunication        = std::shared_ptr<Communication>;
using PtrCommunicationFactory = std::shared_ptr<CommunicationFactory>;
using PtrRequest              = std::shared_ptr<Request>;
} // namespace com
} // namespace precice
