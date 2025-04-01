#pragma once

#include <memory>

namespace precice::mesh {

class Data;
class Group;
class Mesh;
class DataConfiguration;
class MeshConfiguration;

using PtrData              = std::shared_ptr<Data>;
using PtrGroup             = std::shared_ptr<Group>;
using PtrMesh              = std::shared_ptr<Mesh>;
using PtrDataConfiguration = std::shared_ptr<DataConfiguration>;
using PtrMeshConfiguration = std::shared_ptr<MeshConfiguration>;

} // namespace precice::mesh
