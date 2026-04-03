#include "com/IntraCommunication.hpp"
#include "utils/assertion.hpp"

namespace precice::com {

void IntraCommunication::sendRange(precice::span<const double> itemsToSend, Rank rankReceiver)
{
  int size = itemsToSend.size();
  send(size, rankReceiver);
  if (size > 0) {
    send(itemsToSend, rankReceiver);
  }
}

void IntraCommunication::sendRange(precice::span<const int> itemsToSend, Rank rankReceiver)
{
  int size = itemsToSend.size();
  send(size, rankReceiver);
  if (size > 0) {
    send(itemsToSend, rankReceiver);
  }
}

std::vector<int> IntraCommunication::receiveRange(Rank rankSender, AsVectorTag<int>)
{
  int size{-1};
  receive(size, rankSender);
  PRECICE_ASSERT(size >= 0);
  std::vector<int> result;
  if (size > 0) {
    result.resize(size);
    receive(result, rankSender);
  }
  return result;
}

std::vector<double> IntraCommunication::receiveRange(Rank rankSender, AsVectorTag<double>)
{
  int size{-1};
  receive(size, rankSender);
  PRECICE_ASSERT(size >= 0);
  std::vector<double> result;
  if (size > 0) {
    result.resize(size);
    receive(result, rankSender);
  }
  return result;
}

} // namespace precice::com
