#include <Eigen/Core>
#include <algorithm>
#include <limits>
#include <memory>
#include <ostream>
#include <utility>

#include "WatchPoint.hpp"
#include "com/Communication.hpp"
#include "com/SharedPointer.hpp"
#include "logging/LogMacros.hpp"
#include "mapping/Polation.hpp"
#include "mesh/Data.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/Vertex.hpp"
#include "precice/types.hpp"
#include "query/Index.hpp"
#include "utils/MasterSlave.hpp"
#include "utils/assertion.hpp"

namespace precice {
namespace impl {

WatchPoint::WatchPoint(
    Eigen::VectorXd    pointCoords,
    mesh::PtrMesh      meshToWatch,
    const std::string &exportFilename)
    : _point(std::move(pointCoords)),
      _mesh(std::move(meshToWatch)),
      _txtWriter(exportFilename)
{
  PRECICE_ASSERT(_mesh);
  PRECICE_ASSERT(_point.size() == _mesh->getDimensions(), _point.size(),
                 _mesh->getDimensions());

  io::TXTTableWriter::DataType vectorType = _mesh->getDimensions() == 2
                                                ? io::TXTTableWriter::VECTOR2D
                                                : io::TXTTableWriter::VECTOR3D;
  _txtWriter.addData("Time", io::TXTTableWriter::DOUBLE);
  _txtWriter.addData("Coordinate", vectorType);
  for (size_t i = 0; i < _mesh->data().size(); i++) {
    _dataToExport.push_back(_mesh->data()[i]);
    if (_dataToExport[i]->getDimensions() > 1) {
      _txtWriter.addData(_dataToExport[i]->getName(), vectorType);
    } else {
      _txtWriter.addData(_dataToExport[i]->getName(), io::TXTTableWriter::DOUBLE);
    }
  }
}

const mesh::PtrMesh &WatchPoint::mesh() const
{
  return _mesh;
}

void WatchPoint::initialize()
{
  PRECICE_TRACE();

  if (_mesh->vertices().size() > 0) {
    auto match        = query::Index{_mesh}.findNearestProjection(_point, 4);
    _interpolation    = std::make_unique<mapping::Polation>(match.polation);
    _shortestDistance = match.distance;
  }

  if (utils::MasterSlave::isSlave()) {
    utils::MasterSlave::_communication->send(_shortestDistance, 0);
    utils::MasterSlave::_communication->receive(_isClosest, 0);
  }

  if (utils::MasterSlave::isMaster()) {
    int    closestRank           = 0;
    double closestDistanceGlobal = _shortestDistance;
    double closestDistanceLocal  = std::numeric_limits<double>::max();
    for (Rank rankSlave : utils::MasterSlave::allSlaves()) {
      utils::MasterSlave::_communication->receive(closestDistanceLocal, rankSlave);
      if (closestDistanceLocal < closestDistanceGlobal) {
        closestDistanceGlobal = closestDistanceLocal;
        closestRank           = rankSlave;
      }
    }
    _isClosest = closestRank == 0;
    for (Rank rankSlave : utils::MasterSlave::allSlaves()) {
      utils::MasterSlave::_communication->send(closestRank == rankSlave, rankSlave);
    }
  }

  PRECICE_DEBUG("Rank: {}, isClosest: {}", utils::MasterSlave::getRank(), _isClosest);
}

void WatchPoint::exportPointData(
    double time)
{
  if (!_isClosest) {
    return;
  }

  //PRECICE_ASSERT(_vertices.size() == _weights.size());
  _txtWriter.writeData("Time", time);
  // Export watch point coordinates
  Eigen::VectorXd coords = Eigen::VectorXd::Constant(_mesh->getDimensions(), 0.0);
  for (const auto &elem : _interpolation->getWeightedElements()) {
    coords += elem.weight * _mesh->vertices()[elem.vertexID].getCoords();
  }
  if (coords.size() == 2) {
    _txtWriter.writeData("Coordinate", Eigen::Vector2d(coords));
  } else {
    _txtWriter.writeData("Coordinate", Eigen::Vector3d(coords));
  }
  // Export watch point data
  for (auto &elem : _dataToExport) {
    if (elem->getDimensions() > 1) {
      Eigen::VectorXd toExport = Eigen::VectorXd::Zero(_mesh->getDimensions());
      getValue(toExport, elem);
      if (coords.size() == 2) {
        _txtWriter.writeData(elem->getName(), Eigen::Vector2d(toExport));
      } else {
        _txtWriter.writeData(elem->getName(), Eigen::Vector3d(toExport));
      }

    } else {
      double valueToExport = 0.0;
      getValue(valueToExport, elem);
      _txtWriter.writeData(elem->getName(), valueToExport);
    }
  }
}

void WatchPoint::getValue(
    Eigen::VectorXd &value,
    mesh::PtrData &  data)
{
  int                    dim = _mesh->getDimensions();
  Eigen::VectorXd        temp(dim);
  const Eigen::VectorXd &values = data->values();
  for (const auto &elem : _interpolation->getWeightedElements()) {
    int offset = elem.vertexID * dim;
    for (int i = 0; i < dim; i++) {
      temp[i] = values[offset + i];
    }
    temp *= elem.weight;
    value += temp;
  }
}

void WatchPoint::getValue(
    double &       value,
    mesh::PtrData &data)
{
  const Eigen::VectorXd &values = data->values();
  for (const auto &elem : _interpolation->getWeightedElements()) {
    value += elem.weight * values[elem.vertexID];
  }
}

} // namespace impl
} // namespace precice
