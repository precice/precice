#include "partition/Partition.hpp"
#include "mesh/Mesh.hpp"

namespace precice {
namespace partition {

Partition::Partition(mesh::PtrMesh mesh)
    : _mesh(mesh)
{
}

} // namespace partition
} // namespace precice
