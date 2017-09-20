#include "ImportSTL.hpp"
#include "utils/Globals.hpp"
#include "mesh/Mesh.hpp"

namespace precice {
namespace io {

logging::Logger ImportSTL:: _log("io::ImportSTL");

ImportSTL:: ImportSTL
(
  const std::string& location )
:
  Import(location)
{
  // TODO
}


void ImportSTL:: doImport
(
  const std::string& name,
  mesh::Mesh&        mesh )
{
  TRACE(name);
  assertion(mesh.getDimensions() == 3, mesh.getDimensions());
  // TODO
}

}} // namespace precice, io
