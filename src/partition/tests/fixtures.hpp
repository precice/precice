#include "partition/ReceivedPartition.hpp"

namespace precice {
namespace partition {

/*
 * @brief A fixture that is used to access private functions of the receivedPartition class.
 *
 * The fixture can be used to call private functions for individual testing.
 */
struct ReceivedPartitionFixture {
  static void createOwnerInformation(ReceivedPartition &part)
  {
    part.createOwnerInformation();
  }
  static void tagMeshFirstRound(ReceivedPartition &part)
  {
    part.tagMeshFirstRound();
  }
  static void prepareBoundingBox(ReceivedPartition &part)
  {
    part.prepareBoundingBox();
  }
};

} // namespace partition
} // namespace precice
