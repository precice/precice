#include "io/ExportCSV.hpp"

#include <Eigen/Core>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <memory>
#include "io/Export.hpp"
#include "logging/LogMacros.hpp"
#include "mesh/Data.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/Vertex.hpp"
#include "utils/IntraComm.hpp"
#include "utils/assertion.hpp"

namespace precice::io {

namespace {
struct StridedAccess {
  double const *ptr;
  int           stride;

  double operator*() const
  {
    return *ptr;
  }

  void next()
  {
    std::advance(ptr, stride);
  }
};
} // namespace

ExportCSV::ExportCSV(
    std::string_view  participantName,
    std::string_view  location,
    const mesh::Mesh &mesh,
    ExportKind        kind,
    int               frequency,
    int               rank,
    int               size)
    : Export(participantName, location, mesh, kind, frequency, rank, size){};

void ExportCSV::doExport(int index, double time)
{
  PRECICE_TRACE(index, time, _mesh->getName());
  PRECICE_ASSERT(index >= 0);
  PRECICE_ASSERT(time >= 0.0);

  if (!keepExport(index))
    return;

  // Construct filename
  std::string filename;
  if (isParallel()) {
    // Mesh-Participant.r2.it2
    filename = fmt::format("{}-{}.{}_{}.csv", _mesh->getName(), _participantName, _rank, formatIndex(index));
  } else {
    // Mesh-Participant.it2
    filename = fmt::format("{}-{}.{}.csv", _mesh->getName(), _participantName, formatIndex(index));
  }

  namespace fs = std::filesystem;
  fs::path outfile(_location);
  if (not _location.empty()) {
    fs::create_directories(outfile);
  }
  outfile /= filename;

  // Prepare filestream
  std::ofstream outFile(outfile.string(), std::ios::trunc);
  const bool    is3d = (_mesh->getDimensions() == 3);

  // write header
  outFile << "PosX;PosY";
  if (is3d) {
    outFile << ";PosZ";
  }
  outFile << ";Rank";
  for (const auto &data : _mesh->data()) {
    auto dataName = data->getName();
    auto dim      = data->getDimensions();
    outFile << ';' << dataName;
    if (dim == 2) {
      outFile << "X;" << dataName << 'Y';
    } else if (dim == 3) {
      outFile << "X;" << dataName << "Y;" << dataName << 'Z';
    }
  }
  outFile << '\n';

  // Prepare writing data
  std::vector<StridedAccess> dataColumns;
  for (const auto &data : _mesh->data()) {
    if (data->timeStepsStorage().empty()) {
      if (data->hasGradient()) {
        data->timeStepsStorage().setSampleAtTime(0, time::Sample(data->getDimensions(), _mesh->nVertices(), _mesh->getDimensions()).setZero());
      } else {
        data->timeStepsStorage().setSampleAtTime(0, time::Sample(data->getDimensions(), _mesh->nVertices()).setZero());
      }
    }

    auto          dim    = data->getDimensions();
    double const *values = data->timeStepsStorage().last().sample.values.data();
    for (int i = 0; i < dim; ++i) {
      dataColumns.push_back({std::next(values, i), dim});
    }
  }

  // write vertex data
  const std::string rankCol = ";" + std::to_string(_rank);
  const auto        size    = _mesh->nVertices();
  for (std::size_t vid = 0; vid < size; ++vid) {
    const auto &vertex = _mesh->vertex(vid);
    outFile << vertex.coord(0) << ';';
    outFile << vertex.coord(1);
    if (is3d) {
      outFile << ";" << vertex.coord(2);
    }
    outFile << rankCol;
    for (auto &dc : dataColumns) {
      outFile << ';' << *dc;
      dc.next();
    }
    outFile << '\n';
  }
}

void ExportCSV::exportSeries() const
{
  // not supported by paraview
}

} // namespace precice::io
