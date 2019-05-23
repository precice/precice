#pragma once

#include "logging/Logger.hpp"
#include "xml/XMLTag.hpp"
#include "com/SharedPointer.hpp"

namespace precice
{
namespace com
{

/**
 * @brief Configuration for communication channels between server and clients or master and slaves.
 * The communication between two solvers is configured in m2n::M2NConfiguration
 */
class CommunicationConfiguration
{
public:
  virtual ~CommunicationConfiguration() {}

  /// Returns a communication object of given type.
  PtrCommunication createCommunication(const xml::XMLTag &tag) const;

private:
  mutable logging::Logger _log{"com::CommunicationConfiguration"};
};

} // namespace com
} // namespace precice
