#include "partition/Partition.hpp"

#include <utility>

namespace precice {
namespace partition {

Partition::Partition(mesh::PtrMesh mesh)
    : _mesh(std::move(mesh))
{
}

} // namespace partition
} // namespace precice
