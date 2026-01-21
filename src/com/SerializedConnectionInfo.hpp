#pragma once

#include <map>
#include <vector>

#include  "precice/impl/Types.hpp"

namespace precice::com {
class Communication;

namespace serialize {

/// serialized representation of ConnectionInfoMap
class SerializedConnectionInfoMap {
public:
  using ConnectionInfoMap = std::map<Rank, std::string>;

  /** serializes a given ConnectionInfoMap
   *
   * Calls assertValid
   */
  static SerializedConnectionInfoMap serialize(const ConnectionInfoMap &cim);

  /// Builds and returns the connection info map represented by the serialized state
  ConnectionInfoMap toConnectionInfoMap() const;

  /// asserts the content for correctness
  void assertValid() const;

  void send(Communication &communication, int rankReceiver) const;

  /// receives a SerializedConnectionInfoMap and calls assertValid before returning
  static SerializedConnectionInfoMap receive(Communication &communication, int rankSender);

  void broadcastSend(Communication &communication) const;

  /// receives a SerializedConnectionInfoMap and calls assertValid before returning
  static SerializedConnectionInfoMap broadcastReceive(Communication &communication);

private:
  SerializedConnectionInfoMap() = default;

  /// Num entries, Rank0, StringLength0, String0, Rank1, ...
  /// @TODO move to size_t once we changed VertexIDs to size_t
  std::vector<int> content;
};

} // namespace serialize

} // namespace precice::com
