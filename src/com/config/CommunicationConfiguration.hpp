#pragma once

#include <string>
#include "com/SharedPointer.hpp"
#include "logging/Logger.hpp"
#include "xml/XMLTag.hpp"

namespace precice {
namespace xml {
class XMLTag;
} // namespace xml

namespace com {

/**
 * @brief Configuration for communication channels between a master and its slaves.
 * The communication between two solvers is configured in m2n::M2NConfiguration
 */
class CommunicationConfiguration {
public:
  virtual ~CommunicationConfiguration() {}

  /// Returns a communication object of given type.
  PtrCommunication createCommunication(const xml::XMLTag &tag) const;

private:
  mutable logging::Logger _log{"com::CommunicationConfiguration"};
};

} // namespace com
} // namespace precice
