#include "partition/Partition.hpp"

#include <utility>

namespace precice::partition {

Partition::Partition(mesh::PtrMesh mesh)
    : _mesh(std::move(mesh))
{
}

} // namespace precice::partition
