#include "WatchIntegral.hpp"
#include <Eigen/Core>
#include <algorithm>
#include <limits>
#include <memory>
#include <ostream>
#include <utility>
#include "com/Communication.hpp"
#include "logging/LogMacros.hpp"
#include "mesh/Data.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/SharedPointer.hpp"
#include "utils/assertion.hpp"
#include "utils/MasterSlave.hpp"

namespace precice {
namespace impl {

WatchIntegral::WatchIntegral(
    mesh::PtrMesh      meshToWatch,
    const std::string &exportFilename)
    : _mesh(std::move(meshToWatch)),
      _txtWriter(exportFilename)
{
  PRECICE_ASSERT(_mesh);

  if(not utils::MasterSlave::isSlave()){
   _txtWriter.addData("Time", io::TXTTableWriter::DOUBLE);
  }

  for (size_t i = 0; i < _mesh->data().size(); i++) {
    _dataToExport.push_back(_mesh->data()[i]);

    if(not utils::MasterSlave::isSlave()){
      if (_dataToExport[i]->getDimensions() > 1) {

        io::TXTTableWriter::DataType vectorType = _dataToExport[i]->getDimensions() == 2
                                                ? io::TXTTableWriter::VECTOR2D
                                                : io::TXTTableWriter::VECTOR3D;
        _txtWriter.addData(_dataToExport[i]->getName(), vectorType);
      } else {
        _txtWriter.addData(_dataToExport[i]->getName(), io::TXTTableWriter::DOUBLE);
      }
    }
  }

  if(not utils::MasterSlave::isSlave()){
    _txtWriter.addData("SurfaceArea", io::TXTTableWriter::DOUBLE);
  }
}

const mesh::PtrMesh &WatchIntegral::mesh() const
{
  return _mesh;
}

void WatchIntegral::exportIntegralData(
    double time)
{
  if(not utils::MasterSlave::isSlave()){
    _txtWriter.writeData("Time", time);
  }
  
  for(auto &elem : _dataToExport){
    
    int dataDimensions = elem->getDimensions();
    
    if(dataDimensions > 1){
      
      Eigen::VectorXd value = calculateVectorData(elem);
      
      if(utils::MasterSlave::getSize() > 1){
        Eigen::VectorXd valueRecv = Eigen::VectorXd::Zero(dataDimensions);
        utils::MasterSlave::allreduceSum(value.data(), valueRecv.data(), value.size());
        value = valueRecv;
      }

      if(not utils::MasterSlave::isSlave()){
        if(dataDimensions == 2){
          _txtWriter.writeData(elem->getName(), Eigen::Vector2d(value));
        }
        else{
          _txtWriter.writeData(elem->getName(), Eigen::Vector3d(value));
        }
      }
    }
    else{
      double value = calculateScalarData(elem);
      if(utils::MasterSlave::getSize() > 1){
        double valueSum = 0.0;
        utils::MasterSlave::allreduceSum(value, valueSum, 1);
        value = valueSum;
      }
      if(not utils::MasterSlave::isSlave()){
        _txtWriter.writeData(elem->getName(), value);
      }
    }
  }

  double surfaceArea = calculateSurfaceArea();
  if(utils::MasterSlave::getSize() > 1){
    double surfaceAreaSum = 0.0;
    utils::MasterSlave::allreduceSum(surfaceArea, surfaceAreaSum, 1);
    surfaceArea = surfaceAreaSum;
  }

  if(not utils::MasterSlave::isSlave()){
    _txtWriter.writeData("SurfaceArea", surfaceArea);
  }
}

Eigen::VectorXd WatchIntegral::calculateVectorData(mesh::PtrData data){
  
  int dim = data->getDimensions();
  const Eigen::VectorXd& values = data->values();
  Eigen::VectorXd value = Eigen::VectorXd::Zero(dim);

  if(_mesh->edges().empty()){
    for (const auto& vertex : _mesh->vertices()) {
      int offset = vertex.getID() * dim;
      for (int i = 0; i < dim; i++) {
        value[i] += values[offset + i];
      }
    }
  }
  else{ // Connectivity information is given
    // For 2D, connectivity elements are edges
    // Calculate the average and multiply by edge length (midpoint rule)
    if(_mesh->getDimensions() == 2){   
      for(const auto& edge : _mesh->edges()){
          
        const auto& vertex1 = edge.vertex(0);
        const auto& vertex2 = edge.vertex(1);
          
        for (int i = 0; i < dim; i++) {
          value[i] += 0.5 * (values[vertex1.getID() * dim + i] + values[vertex2.getID() * dim + i]) * edge.getLength();
        }
      }
    }
    else{ // For 3D, connectivity elements are faces, calculate the average and multply by face area
      for(const auto& face : _mesh->triangles()){
        for(int i = 0; i < dim; ++i){
          value[i] += face.getArea() * (values[face.vertex(0).getID()*dim + i] + values[face.vertex(1).getID()*dim + i] + values[face.vertex(2).getID()*dim + i]) / 3.0;
        }
      }
    }
  }
  return std::move(value);
}

double WatchIntegral::calculateScalarData(mesh::PtrData data){
  
  const Eigen::VectorXd& values = data->values();
  double value = 0.0;

  if(_mesh->edges().empty()){
    for(const auto& vertex : _mesh->vertices()){
      value += values[vertex.getID()];
    }
  }
  else{ // Connectivity information is given
    // For 2D, connectivity elements are edges
    // Calculate the average and multiply by edge length (midpoint rule)
    if(_mesh->getDimensions() == 2){   
      for(const auto& edge : _mesh->edges()){
          
        const auto& vertex1 = edge.vertex(0);
        const auto& vertex2 = edge.vertex(1);
          
        value += 0.5 * (values[vertex1.getID()] + values[vertex2.getID()]) * edge.getLength();
      }
    }
    else{ // For 3D, connectivity elements are faces, calculate the average and multply by face area
      for(const auto& face : _mesh->triangles()){
        double localValue = values[face.vertex(0).getID()] + values[face.vertex(1).getID()] + values[face.vertex(2).getID()];
        value += (localValue * face.getArea()) / 3.0;
      }
    }
  }
  return value;
}

double WatchIntegral::calculateSurfaceArea(){
  
  double surfaceArea = 0.0;

  // 3D and has face connectivity
  if(not _mesh->triangles().empty()){
    for(const auto& face : _mesh->triangles()){
      surfaceArea += face.getArea();
    }
  }
  // 2D and has edge connectivity
  else if(not _mesh->edges().empty()){
    for(const auto& edge : _mesh->edges()){
      surfaceArea += edge.getLength();
    }
  }
  // No connectivity, dimension is not important
  else{
    surfaceArea += static_cast<double>(_mesh->vertices().size());
  }
  return surfaceArea;
}

} // namespace impl
} // namespace precice
