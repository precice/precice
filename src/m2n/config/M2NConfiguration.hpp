#pragma once

#include "logging/Logger.hpp"
#include "m2n/SharedPointer.hpp"
#include "xml/XMLTag.hpp"

#include <tuple>

#include <memory>
#include <string>
#include <vector>

namespace precice
{
namespace m2n
{

/// Configuration for communication channels between solvers.
class M2NConfiguration : public xml::XMLTag::Listener
{
public:
  using SharedPointer = std::shared_ptr<M2NConfiguration>;
  using M2NTuple      = std::tuple<m2n::PtrM2N, std::string, std::string>;

public:
  explicit M2NConfiguration(xml::XMLTag &parent);

  virtual ~M2NConfiguration() {}

  /**
    * @brief Returns the communication object for the given user names.
    *
    * Exits with an error message, when no object is configured for the given
    * user names.
    */
  m2n::PtrM2N getM2N(
      const std::string &from,
      const std::string &to);

  /// Returns all configured communication objects.
  std::vector<M2NTuple> &m2ns()
  {
    return _m2ns;
  }

  virtual void xmlTagCallback(const xml::ConfigurationContext& context, xml::XMLTag &callingTag);

  virtual void xmlEndTagCallback(const xml::ConfigurationContext& context, xml::XMLTag &callingTag) {}

private:
  logging::Logger _log{"m2n::M2NConfiguration"};

  const std::string TAG                     = "m2n";
  const std::string ATTR_DISTRIBUTION_TYPE  = "distribution-type";
  const std::string ATTR_EXCHANGE_DIRECTORY = "exchange-directory";

  const std::string VALUE_GATHER_SCATTER = "gather-scatter";
  const std::string VALUE_POINT_TO_POINT = "point-to-point";

  std::vector<M2NTuple> _m2ns;

  void checkDuplicates(
      const std::string &from,
      const std::string &to);
};

} // namespace m2n
} // namespace precice
