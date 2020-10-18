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
      
      Eigen::VectorXd value = Eigen::VectorXd::Zero(dataDimensions);
      double surfaceArea = 0.0;
      calculateVectorData(value, surfaceArea, elem);
      
      if(utils::MasterSlave::getSize() > 1){
        Eigen::VectorXd valueRecv = Eigen::VectorXd::Zero(dataDimensions);
        utils::MasterSlave::allreduceSum(value.data(), valueRecv.data(), value.size());
        value = valueRecv;
        utils::MasterSlave::allreduceSum(surfaceArea, surfaceArea, 1);
      }

      if(not utils::MasterSlave::isSlave()){
        if(dataDimensions == 2){
          _txtWriter.writeData(elem->getName(), Eigen::Vector2d(value));
        }
        else{
          _txtWriter.writeData(elem->getName(), Eigen::Vector3d(value));
        }
        _txtWriter.writeData("SurfaceArea", surfaceArea);
      }
    
    }
    else{

      double value = 0.0;
      double surfaceArea = 0.0;
      calculateScalarData(value, surfaceArea, elem);

      if(utils::MasterSlave::getSize() > 1){
        utils::MasterSlave::allreduceSum(value, value, 1);
        utils::MasterSlave::allreduceSum(surfaceArea, surfaceArea, 1);
      }
      if(not utils::MasterSlave::isSlave()){
        _txtWriter.writeData(elem->getName(), value);
        _txtWriter.writeData("SurfaceArea", surfaceArea);
      }
    }
  }
}

void WatchIntegral::calculateVectorData(Eigen::VectorXd& value, double& surfaceArea, mesh::PtrData data){
  
  int dim = data->getDimensions();
  const Eigen::VectorXd& values = data->values();
  Eigen::VectorXd temp(dim);
  value = Eigen::VectorXd::Zero(data->getDimensions());
  surfaceArea = 0.0;

  if(_mesh->edges().empty()){
    
    int index = 0;

    for (const auto& vertex : _mesh->vertices()) {
      int offset = vertex.getID() * dim;
      for (int i = 0; i < dim; i++) {
        temp[i] = values[offset + i];
      }
      value += temp;
      index++;
    }
    surfaceArea = static_cast<double>(_mesh->vertices().size());
  }
  else{ // Connectivity information is given
    // For 2D, connectivity elements are edges
    // Calculate the average and multiply by edge length (midpoint rule)
    if(_mesh->getDimensions() == 2){   
      for(const auto& edge : _mesh->edges()){
          
        const auto& vertex1 = edge.vertex(0);
        const auto& vertex2 = edge.vertex(1);
          
        for (int i = 0; i < dim; i++) {
          temp[i] = 0.5 * (values[vertex1.getID() * dim + i] + values[vertex2.getID() * dim + i]) * edge.getLength();
        }

        value += temp;
        surfaceArea += edge.getLength();
      }
    }
    else{ // For 3D, connectivity elements are faces, calculate the average and multply by face area
      for(const auto& face : _mesh->triangles()){
        
        for(int i = 0; i < dim; ++i){
          temp[i] = values[face.vertex(0).getID()*dim + i] + values[face.vertex(1).getID()*dim + i] + values[face.vertex(2).getID()*dim + i];
        }
        
        value += (temp * face.getArea()) / 3.0;
        surfaceArea += face.getArea();
      }
    }
  }
}

void WatchIntegral::calculateScalarData(double& value, double& surfaceArea, mesh::PtrData data){
  
  const Eigen::VectorXd& values = data->values();
  value = 0.0;
  surfaceArea = 0.0;

  if(_mesh->edges().empty()){
    for(const auto& vertex : _mesh->vertices()){
      value += values[vertex.getID()];
    }
    surfaceArea = static_cast<double>(_mesh->vertices().size());
  }
  else{ // Connectivity information is given
    // For 2D, connectivity elements are edges
    // Calculate the average and multiply by edge length (midpoint rule)
    if(_mesh->getDimensions() == 2){   
      for(const auto& edge : _mesh->edges()){
          
        const auto& vertex1 = edge.vertex(0);
        const auto& vertex2 = edge.vertex(1);
          
        value += 0.5 * (values[vertex1.getID()] + values[vertex2.getID()]) * edge.getLength();
        surfaceArea += edge.getLength();
      }
    }
    else{ // For 3D, connectivity elements are faces, calculate the average and multply by face area
      for(const auto& face : _mesh->triangles()){
        double localValue = values[face.vertex(0).getID()] + values[face.vertex(1).getID()] + values[face.vertex(2).getID()];
        value += (localValue * face.getArea()) / 3.0;
        surfaceArea += face.getArea();
      }
    }
  }
}

} // namespace impl
} // namespace precice
