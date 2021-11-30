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
#include "utils/MasterSlave.hpp"
#include "utils/assertion.hpp"

namespace precice {
namespace io {

void ExportCSV::doExport(
    const std::string &name,
    const std::string &location,
    const mesh::Mesh & mesh)
{
  PRECICE_TRACE(name, location, mesh.getName());
  PRECICE_ASSERT(!name.empty());

  // Ignore empty meshes
  if (mesh.vertices().empty()) {
    return;
  }

  // Construct full filename
  std::string filename{name};
  int         rank{0};
  if (utils::MasterSlave::isParallel()) {
    rank = utils::MasterSlave::getRank();
    filename.append("_").append(std::to_string(rank));
  }
  filename.append(".csv");

  namespace fs = boost::filesystem;
  fs::path outfile(location);
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
    auto name = data->getName();
    auto dim  = data->getDimensions();
    outFile << ';' << name;
    if (dim == 2) {
      outFile << "X;" << name << 'Y';
    } else if (dim == 3) {
      outFile << "X;" << name << "Y;" << name << 'Z';
    }
  }
  outFile << '\n';

  // write vertex data
  const std::string rankCol = ";" + std::to_string(rank);
  const auto        size    = mesh.vertices().size();
  for (std::size_t vid = 0; vid < size; ++vid) {
    const auto &vertex = mesh.vertices()[vid];
    outFile << vertex.getCoords()[0] << ';';
    outFile << vertex.getCoords()[1];
    if (is3d) {
      outFile << ";" << vertex.getCoords()[2];
    }
    outFile << rankCol;
    for (const auto &data : mesh.data()) {
      auto dim    = data->getDimensions();
      auto offset = vid * dim;

      outFile << ';' << data->values()[offset];
      if (dim == 1) {
        continue;
      }
      outFile << ';' << data->values()[offset + 1];
      if (dim == 2) {
        continue;
      }
      outFile << ';' << data->values()[offset + 2];
    }
    outFile << '\n';
  }
}

} // namespace io
} // namespace precice
