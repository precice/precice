#include "partition/ReceivedPartition.hpp"

namespace precice {
namespace partition {

/*
 * @brief A fixture that is used to access private functions of the receivedPartition class.
 *
 * The fixture can be used to call private functions for individual testing. 
 */
struct receivedPartitionFixture {
  void createOwnerInformation(ReceivedPartition &part)
  {
    part.createOwnerInformation();
  }
  void tagMeshFirstRound(ReceivedPartition &part)
  {
    part.tagMeshFirstRound();
  }
  void prepareBoundingBox(ReceivedPartition &part)
  {
    part.prepareBoundingBox();
  }
};

} // namespace partition
} // namespace precice
