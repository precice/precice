#include "ImportSTL.hpp"
#include "utils/Globals.hpp"
#include "mesh/Mesh.hpp"

namespace precice {
namespace io {

tarch::logging::Log ImportSTL:: _log("precice::io::ImportSTL");

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
  preciceTrace1("doImport()", name);
  assertion1(mesh.getDimensions() == 3, mesh.getDimensions());
  // TODO
}

}} // namespace precice, io
