#include "WatchIntegral.hpp"
#include <Eigen/Core>
#include "logging/LogMacros.hpp"
#include "mesh/Data.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/Vertex.hpp"
#include "utils/MasterSlave.hpp"

namespace precice {
namespace impl {

WatchIntegral::WatchIntegral(
    mesh::PtrMesh      meshToWatch,
    const std::string &exportFilename,
    bool               isScalingOn)
    : _mesh(std::move(meshToWatch)),
      _txtWriter(exportFilename),
      _isScalingOn(isScalingOn)
{
  PRECICE_ASSERT(_mesh);

  if (not utils::MasterSlave::isSlave()) {
    _txtWriter.addData("Time", io::TXTTableWriter::DOUBLE);
  }

  for (size_t i = 0; i < _mesh->data().size(); i++) {
    _dataToExport.push_back(_mesh->data()[i]);

    if (not utils::MasterSlave::isSlave()) {
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
}

void WatchIntegral::initialize()
{
  // Do not add surface area column if there is no connectivity
  if ((not utils::MasterSlave::isSlave()) and (not _mesh->edges().empty())) {
    _txtWriter.addData("SurfaceArea", io::TXTTableWriter::DOUBLE);
  }
  if (_isScalingOn and (_mesh->edges().empty())) {
    PRECICE_WARN("Watch-integral is configured with scaling option on; however, mesh " << _mesh->getName() << " does not contain connectivity information. Therefore, the integral will be calculated without scaling.");
  }
}

void WatchIntegral::exportIntegralData(
    double time)
{

  if (not utils::MasterSlave::isSlave()) {
    _txtWriter.writeData("Time", time);
  }

  for (auto &elem : _dataToExport) {

    int dataDimensions = elem->getDimensions();

    if (dataDimensions > 1) {

      Eigen::VectorXd value = calculateVectorData(elem);

      if (utils::MasterSlave::getSize() > 1) {
        Eigen::VectorXd valueRecv = Eigen::VectorXd::Zero(dataDimensions);
        utils::MasterSlave::reduceSum(value.data(), valueRecv.data(), value.size());
        value = valueRecv;
      }

      if (not utils::MasterSlave::isSlave()) {
        if (dataDimensions == 2) {
          _txtWriter.writeData(elem->getName(), Eigen::Vector2d(value));
        } else {
          _txtWriter.writeData(elem->getName(), Eigen::Vector3d(value));
        }
      }
    } else {
      double value = calculateScalarData(elem);
      if (utils::MasterSlave::getSize() > 1) {
        double valueSum = 0.0;
        utils::MasterSlave::reduceSum(&value, &valueSum, 1);
        value = valueSum;
      }
      if (not utils::MasterSlave::isSlave()) {
        _txtWriter.writeData(elem->getName(), value);
      }
    }
  }

  // Calculate surface area only if there is connectivity information
  if (not _mesh->edges().empty()) {
    double surfaceArea = calculateSurfaceArea();
    if (utils::MasterSlave::getSize() > 1) {
      double surfaceAreaSum = 0.0;
      utils::MasterSlave::reduceSum(&surfaceArea, &surfaceAreaSum, 1);
      surfaceArea = surfaceAreaSum;
    }
    if (not utils::MasterSlave::isSlave()) {
      _txtWriter.writeData("SurfaceArea", surfaceArea);
    }
  } else {
    // Empty partitions should also call reduceSum to prevent deadlock
    if (utils::MasterSlave::getSize() > 1) {
      double surfaceArea    = 0.0;
      double surfaceAreaSum = 0.0;
      utils::MasterSlave::reduceSum(&surfaceArea, &surfaceAreaSum, 1);
    }
  }
}

Eigen::VectorXd WatchIntegral::calculateVectorData(mesh::PtrData data)
{

  int                    dim    = data->getDimensions();
  const Eigen::VectorXd &values = data->values();
  Eigen::VectorXd        sum    = Eigen::VectorXd::Zero(dim);

  if (_mesh->edges().empty() || (not _isScalingOn)) {
    for (const auto &vertex : _mesh->vertices()) {
      int offset = vertex.getID() * dim;
      for (int i = 0; i < dim; i++) {
        sum[i] += values[offset + i];
      }
    }
  } else { // Connectivity information is given
    // For 2D, connectivity elements are edges
    // Calculate the average and multiply by edge length (midpoint rule)
    if (_mesh->getDimensions() == 2) {
      for (const auto &edge : _mesh->edges()) {

        const auto &vertex1 = edge.vertex(0);
        const auto &vertex2 = edge.vertex(1);

        for (int i = 0; i < dim; i++) {
          sum[i] += 0.5 * (values[vertex1.getID() * dim + i] + values[vertex2.getID() * dim + i]) * edge.getLength();
        }
      }
    } else { // For 3D, connectivity elements are faces, calculate the average and multply by face area
      for (const auto &face : _mesh->triangles()) {
        for (int i = 0; i < dim; ++i) {
          sum[i] += face.getArea() * (values[face.vertex(0).getID() * dim + i] + values[face.vertex(1).getID() * dim + i] + values[face.vertex(2).getID() * dim + i]) / 3.0;
        }
      }
    }
  }
  return sum;
}

double WatchIntegral::calculateScalarData(mesh::PtrData data)
{

  const Eigen::VectorXd &values = data->values();
  double                 sum    = 0.0;

  if (_mesh->edges().empty() || not _isScalingOn) {
    for (const auto &vertex : _mesh->vertices()) {
      sum += values[vertex.getID()];
    }
  } else { // Connectivity information is given
    // For 2D, connectivity elements are edges
    // Calculate the average and multiply by edge length (midpoint rule)
    if (_mesh->getDimensions() == 2) {
      for (const auto &edge : _mesh->edges()) {

        const auto &vertex1 = edge.vertex(0);
        const auto &vertex2 = edge.vertex(1);

        sum += 0.5 * (values[vertex1.getID()] + values[vertex2.getID()]) * edge.getLength();
      }
    } else { // For 3D, connectivity elements are faces, calculate the average and multply by face area
      for (const auto &face : _mesh->triangles()) {
        double localValue = values[face.vertex(0).getID()] + values[face.vertex(1).getID()] + values[face.vertex(2).getID()];
        sum += (localValue * face.getArea()) / 3.0;
      }
    }
  }
  return sum;
}

double WatchIntegral::calculateSurfaceArea()
{

  PRECICE_ASSERT(not _mesh->edges().empty());
  double surfaceArea = 0.0;
  if (_mesh->getDimensions() == 3) {
    for (const auto &face : _mesh->triangles()) {
      surfaceArea += face.getArea();
    }
  } else {
    for (const auto &edge : _mesh->edges()) {
      surfaceArea += edge.getLength();
    }
  }
  return surfaceArea;
}

} // namespace impl
} // namespace precice
