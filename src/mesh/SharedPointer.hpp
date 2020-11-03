#pragma once

#include <memory>

namespace precice {
namespace mesh {

class Data;
class Gradient;
class Group;
class Mesh;
class DataConfiguration;
class GradientConfiguration;
class MeshConfiguration;

using PtrData                  = std::shared_ptr<Data>;
using PtrGradient              = std::shared_ptr<Gradient>;
using PtrGroup                 = std::shared_ptr<Group>;
using PtrMesh                  = std::shared_ptr<Mesh>;
using PtrDataConfiguration     = std::shared_ptr<DataConfiguration>;
using PtrGradientConfiguration = std::shared_ptr<GradientConfiguration>;
using PtrMeshConfiguration     = std::shared_ptr<MeshConfiguration>;

} // namespace mesh
} // namespace precice
