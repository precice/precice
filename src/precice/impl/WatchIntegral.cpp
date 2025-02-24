#include "WatchIntegral.hpp"
#include <Eigen/Core>
#include "logging/LogMacros.hpp"
#include "mesh/Data.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/Utils.hpp"
#include "mesh/Vertex.hpp"
#include "utils/IntraComm.hpp"

namespace precice::impl {

WatchIntegral::WatchIntegral(
    mesh::PtrMesh      meshToWatch,
    const std::string &exportFilename,
    bool               isScalingOn)
    : _mesh(std::move(meshToWatch)),
      _txtWriter(exportFilename),
      _isScalingOn(isScalingOn)
{
  PRECICE_ASSERT(_mesh);

  if (not utils::IntraComm::isSecondary()) {
    _txtWriter.addData("Time", io::TXTTableWriter::DOUBLE);
  }

  for (size_t i = 0; i < _mesh->data().size(); i++) {
    _dataToExport.push_back(_mesh->data()[i]);

    if (not utils::IntraComm::isSecondary()) {
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
  if ((not utils::IntraComm::isSecondary()) and (not _mesh->edges().empty())) {
    _txtWriter.addData("SurfaceArea", io::TXTTableWriter::DOUBLE);
  }
  PRECICE_WARN_IF(_isScalingOn && _mesh->edges().empty(),
                  "Watch-integral is configured with scaling option on; "
                  "however, mesh {} does not contain connectivity information. "
                  "Therefore, the integral will be calculated without scaling.",
                  _mesh->getName());
}

void WatchIntegral::exportIntegralData(
    double time)
{

  if (not utils::IntraComm::isSecondary()) {
    _txtWriter.writeData("Time", time);
  }

  for (auto &elem : _dataToExport) {
    const int dataDimensions = elem->getDimensions();
    auto      integral       = calculateIntegral(time, elem);

    if (utils::IntraComm::getSize() > 1) {
      Eigen::VectorXd valueRecv = Eigen::VectorXd::Zero(dataDimensions);
      utils::IntraComm::reduceSum(integral, valueRecv);
      integral = std::move(valueRecv);
    }
    if (not utils::IntraComm::isSecondary()) {
      if (dataDimensions == 1) {
        _txtWriter.writeData(elem->getName(), integral[0]);
      } else if (dataDimensions == 2) {
        _txtWriter.writeData(elem->getName(), Eigen::Vector2d(integral));
      } else {
        _txtWriter.writeData(elem->getName(), Eigen::Vector3d(integral));
      }
    }
  }

  // Calculate surface area only if there is connectivity information
  if (not _mesh->edges().empty()) {
    double surfaceArea = calculateSurfaceArea();
    if (utils::IntraComm::getSize() > 1) {
      double surfaceAreaSum = 0.0;
      utils::IntraComm::reduceSum(surfaceArea, surfaceAreaSum);
      surfaceArea = surfaceAreaSum;
    }
    if (not utils::IntraComm::isSecondary()) {
      _txtWriter.writeData("SurfaceArea", surfaceArea);
    }
  } else {
    // Empty partitions should also call reduceSum to prevent deadlock
    if (utils::IntraComm::getSize() > 1) {
      double surfaceArea    = 0.0;
      double surfaceAreaSum = 0.0;
      utils::IntraComm::reduceSum(surfaceArea, surfaceAreaSum);
    }
  }
}

Eigen::VectorXd WatchIntegral::calculateIntegral(double time, const mesh::PtrData &data) const
{
  int             dim    = data->getDimensions();
  const auto      sample = data->timeStepsStorage().sample(time);
  Eigen::VectorXd sum    = Eigen::VectorXd::Zero(dim);

  if (_mesh->edges().empty() || (not _isScalingOn)) {
    for (const auto &vertex : _mesh->vertices()) {
      int offset = vertex.getID() * dim;
      for (int i = 0; i < dim; i++) {
        sum[i] += sample(offset + i);
      }
    }
    return sum;
  } else { // Connectivity information is given
    return mesh::integrateSurface(_mesh, data->timeStepsStorage().sample(time).values());
  }
}

double WatchIntegral::calculateSurfaceArea() const
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

} // namespace precice::impl
