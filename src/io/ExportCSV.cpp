#include "io/ExportCSV.hpp"

#include <Eigen/Core>
#include <boost/filesystem.hpp>
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

void ExportCSV::doExport(
    const std::string &name,
    const std::string &location,
    const mesh::Mesh & mesh)
{
  PRECICE_TRACE(name, location, mesh.getName());
  PRECICE_ASSERT(!name.empty());

  // Ignore empty meshes
  if (mesh.empty()) {
    return;
  }

  // Construct full filename
  std::string filename{name};
  int         rank{0};
  if (utils::IntraComm::isParallel()) {
    rank = utils::IntraComm::getRank();
    filename.append("_").append(std::to_string(rank));
  }
  filename.append(".csv");

  namespace fs = boost::filesystem;
  fs::path outfile(location);
  if (not location.empty()) {
    fs::create_directories(outfile);
  }
  outfile /= filename;

  // Prepare filestream
  std::ofstream outFile(outfile.string(), std::ios::trunc);
  const bool    is3d = (mesh.getDimensions() == 3);

  // write header
  outFile << "PosX;PosY";
  if (is3d) {
    outFile << ";PosZ";
  }
  outFile << ";Rank";
  for (const auto &data : mesh.data()) {
    auto dataName = data->getName();
    auto dim      = data->getDimensions();
    PRECICE_ASSERT(static_cast<std::size_t>(data->values().size()) == mesh.nVertices() * dim);
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
  for (const auto &data : mesh.data()) {
    auto    dim    = data->getDimensions();
    double *values = data->values().data();
    for (int i = 0; i < dim; ++i) {
      dataColumns.push_back({std::next(values, i), dim});
    }
  }

  // write vertex data
  const std::string rankCol = ";" + std::to_string(rank);
  const auto        size    = mesh.nVertices();
  for (std::size_t vid = 0; vid < size; ++vid) {
    const auto &vertex = mesh.vertices()[vid];
    outFile << vertex.getCoords()[0] << ';';
    outFile << vertex.getCoords()[1];
    if (is3d) {
      outFile << ";" << vertex.getCoords()[2];
    }
    outFile << rankCol;
    for (auto &dc : dataColumns) {
      outFile << ';' << *dc;
      dc.next();
    }
    outFile << '\n';
  }
}

} // namespace precice::io
