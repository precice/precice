#pragma once

#include <numeric>

#include "com/SerializedConnectionInfo.hpp"
#include "utils/assertion.hpp"

namespace precice::com::serialize {

// SerializedConnectionInfoMap

SerializedConnectionInfoMap SerializedConnectionInfoMap::serialize(const ConnectionInfoMap &cim)
{
  SerializedConnectionInfoMap scim;
  auto                       &content = scim.content;

  auto elements = std::transform_reduce(
      cim.begin(), cim.end(),
      1, // num entries
      std::plus<size_t>{}, [](const auto &kv) {
        // Rank, size, characters
        return 2 + kv.second.size();
      });
  content.reserve(elements);

  content.push_back(cim.size());
  for (const auto &[rank, string] : cim) {
    content.push_back(rank);
    content.push_back(string.size());
    content.insert(content.end(), string.begin(), string.end());
  }
  PRECICE_ASSERT(content.size() == static_cast<std::size_t>(elements));
  scim.assertValid();
  return scim;
}

SerializedConnectionInfoMap::ConnectionInfoMap SerializedConnectionInfoMap::toConnectionInfoMap() const
{
  auto begin      = content.begin();
  int  numEntries = *begin;

  if (numEntries == 0) {
    return {};
  }
  std::advance(begin, 1);

  ConnectionInfoMap cim;
  for (int entry = 0; entry < numEntries; ++entry) {
    // offset 0: rank
    auto rank = *begin;
    std::advance(begin, 1);
    // offset 1: string size n
    auto size = *begin;
    std::advance(begin, 1);
    // offset 2: n characters
    auto stringStart = begin;
    std::advance(begin, size);

    cim.emplace(rank, std::string(stringStart, begin));
  }
  return cim;
}

void SerializedConnectionInfoMap::assertValid() const
{
  // TODO
}

void SerializedConnectionInfoMap::send(Communication &communication, int rankReceiver) const
{
  communication.sendRange(content, rankReceiver);
}

SerializedConnectionInfoMap SerializedConnectionInfoMap::receive(Communication &communication, int rankSender)
{
  SerializedConnectionInfoMap scim;
  scim.content = communication.receiveRange(rankSender, asVector<int>);
  scim.assertValid();
  return scim;
}

void SerializedConnectionInfoMap::broadcastSend(Communication &communication) const
{
  communication.broadcast(content);
}

SerializedConnectionInfoMap SerializedConnectionInfoMap::broadcastReceive(Communication &communication)
{
  SerializedConnectionInfoMap scim;
  communication.broadcast(scim.content, 0);
  scim.assertValid();
  return scim;
}

} // namespace precice::com::serialize
