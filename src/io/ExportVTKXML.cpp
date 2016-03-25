#include "ExportVTKXML.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/Quad.hpp"
#include "utils/Globals.hpp"
#include "utils/String.hpp"
#include "utils/Parallel.hpp"
#include "tarch/plotter/griddata/unstructured/vtk/VTKTextFileWriter.h"
#include "tarch/la/WrappedVector.h"
#include "Eigen/Dense"
#include <iostream>
#include <string>
#include <fstream>

namespace precice {
namespace io {

tarch::logging::Log ExportVTKXML:: _log("precice::io::ExportVTKXML");

ExportVTKXML:: ExportVTKXML
(
  bool writeNormals,
  bool parallelWrite  )
:
  Export(),
  _writeNormals(writeNormals),
  _parallelWrite(parallelWrite)
{}

//int ExportVTKXML:: getType() const
//{
//  return constants::exportVTK();
//}

void ExportVTKXML:: doExport
(
  const std::string& filename,
  mesh::Mesh&        mesh)
{
	std::ofstream outFile;
	if (_parallelWrite) { // if parallel write enabled
		if (utils::Parallel::getCommunicatorWorld() != -1) { // if precice in parallel mode
			// code for parallel write for parallel precice
			int rank = utils::Parallel::getLocalProcessRank();
			if (rank == 0) {
				writeMasterFile(filename + "_master.pvtu", mesh);
			}
		} else {
			// code for parallel write of serial precice
		}
	} else {
		if (utils::Parallel::getCommunicatorWorld() != -1) { // if precice in parallel mode
			// code for serial write for parallel precice
			// maybe throw error saying inefficient?
		} else {
			// code for serial write of serial precice
		}
	}
}

void ExportVTKXML::writeMasterFile
(
  const std::string& filename,
  mesh::Mesh&        mesh)
{
  std::ofstream outMasterFile;
  std::string fullMasterFilename(filename);
  //utils::checkAppendExtension(fullMasterFilename, ".pvtu");

  outMasterFile.open(fullMasterFilename.c_str());
  preciceCheck(outMasterFile, "doExport()", "Could not open master file \"" << fullMasterFilename
      << "\" for VTKXML export!");

  outMasterFile << "<?xml version=\"1.0\"?>" << std::endl;
  outMasterFile << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"";
  outMasterFile << (utils::isMachineBigEndian() ? "\"BigEndian\">" : "\"LittleEndian\">")  << std::endl;
  outMasterFile << "   <PUnstructuredGrid GhostLevel=\"0\">" << std::endl;

  outMasterFile << "      <PPoints>" << std::endl;
  outMasterFile << "         <PDataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"" << mesh.getDimensions() << "\"/>" << std::endl;
  outMasterFile << "      </PPoints>" << std::endl;

  outMasterFile << "      <PCells>" << std::endl;
  outMasterFile << "         <PDataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\"/>" << std::endl;
  outMasterFile << "         <PDataArray type=\"Int32\" Name=\"offsets\"      NumberOfComponents=\"1\"/>" << std::endl;
  outMasterFile << "         <PDataArray type=\"UInt8\" Name=\"types\"        NumberOfComponents=\"1\"/>" << std::endl;
  outMasterFile << "      </PCells>" << std::endl;

  std::vector<std::string> scalarDataNames("");
  std::vector<std::string> vectorDataNames("");
  std::vector<int> vectorDataDimensions;
  if (_writeNormals) {
	  vectorDataNames.push_back("VertexNormals ");
	  vectorDataDimensions.push_back(3);
  }
  for (mesh::PtrData data : mesh.data()) {
	  int dimensions = data->getDimensions();
	  std::string dataName = data->getName();
	  if ( dimensions == 1) {
		  scalarDataNames.push_back(dataName);
	  } else if ( dimensions > 1 ) {
		  vectorDataNames.push_back(dataName + " ");
		  vectorDataDimensions.push_back(dimensions);
	  } else {
		  // TODO: Change to an assertion later
		  std::cout << "Error! Dimensions <= 0";
		  exit(-1);
	  }
  }

  // write scalar data names
  outMasterFile << "      <PPointData Scalars=\"";
  for (std::vector<std::string>::iterator name = scalarDataNames.begin(); name != scalarDataNames.end(); ++name) {
	  outMasterFile << name << " ";
  }
  // write vector data names
  outMasterFile << "\" Vectors=\"";
  for (std::vector<std::string>::iterator name = vectorDataNames.begin(); name != vectorDataNames.end(); ++name) {
	  outMasterFile << name << " ";
  }
  outMasterFile << "\">" << std::endl;

  for (std::vector<std::string>::iterator name = scalarDataNames.begin(); name != scalarDataNames.end(); ++name) {
	  outMasterFile << "         <PDataArray type=\"Int32\" Name=\""<< name << "\" NumberOfComponents=\"" << 1 << "\"/>" << std::endl;
  }
  int i = 0;
  for (std::vector<std::string>::iterator name = vectorDataNames.begin(); name != vectorDataNames.end(); ++name) {
	  outMasterFile << "         <PDataArray type=\"Int32\" Name=\""<< name << "\" NumberOfComponents=\"" << vectorDataDimensions[i] << "\"/>" << std::endl;
	  ++i;
  }
  outMasterFile << "      </PPointData>" << std::endl;

  i = 0;
  for (; i < utils::Parallel::getCommunicatorSize(); i++) {
	  outMasterFile << "      <Piece Source=\"output_proc" << i << ".vtu\"/>" << std::endl;
  }

  outMasterFile << "   </PUnstructuredGrid>" << std::endl;
  outMasterFile << "</VTKFile>" << std::endl;

  outMasterFile.close();
}

void ExportVTKXML::exportGeometry
(
  std::ofstream& outFile,
  mesh::Mesh&    mesh)
{

}

void ExportVTKXML:: exportData
(
  std::ofstream& outFile,
  mesh::Mesh&    mesh)
{

}

void ExportVTKXML:: initializeWriting
(
  const std::string& filename,
  std::ofstream&     filestream)
{

}

void ExportVTKXML:: writeHeader
(
  std::ostream& outFile)
{

}

void ExportVTKXML:: writeVertex
(
  const utils::DynVector& position,
  std::ostream&           outFile)
{

}


void ExportVTKXML:: writeTriangle
(
  int           vertexIndices[3],
  std::ostream& outFile)
{

}

void ExportVTKXML:: writeQuadrangle
(
  int           vertexIndices[4],
  std::ostream& outFile)
{

}

void ExportVTKXML:: writeLine
(
  int           vertexIndices[2],
  std::ostream& outFile)
{

}

}} // namespace precice, io

