#pragma once
#include "mesh/Mesh.hpp"

namespace precice{
namespace mesh{

class BoundingBox{

    public:
        BoundingBox(const Mesh& );
        void prepareBoundingBox();
        bool isVertexInBB(const mesh::Vertex &vertex);
    private:
        const Mesh& _mesh;
        const int _dimensions;
        std::vector<std::pair<double, double>> _bounds;

};

} // namespace mesh
} // namespace precice

