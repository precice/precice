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
  _parallelWrite(parallelWrite),
  _isDataNamesAndDimensionsProcessed(false),
  _isCellPresent(),
  _meshDimensions(3)
{
	// TODO: Not working: find a way to find if precice is in parallel mode.
#ifndef PRECICE_NO_MPI
  _isProgramParallel = false;
#else
  _isProgramParallel = true;
#endif
}

int ExportVTKXML:: getType() const
{
  return constants::exportVTKXML();
}

void ExportVTKXML:: doExport
(
  const std::string& filename,
  mesh::Mesh&        mesh)
{
  std::ofstream outFile;
  //preciceDebug("_parallelWrite: " + std::to_string(_parallelWrite));
  //preciceDebug("_isProgramParallel: " + std::to_string(_isProgramParallel));

  if (_parallelWrite) { // if parallel write enabled
    if (1) { // if precice in parallel mode
      // code for parallel write for parallel precice
      //preciceDebug("Writing Parallel File..");
      processDataNamesAndDimensions(mesh);
      int rank = utils::Parallel::getProcessRank();
      if (rank == 0) {
        writeMasterFile(filename, mesh);
      }
      writeSubFile(filename, mesh);
    } else {
      // code for parallel write of serial precice
    }
  } else {
    if (_isProgramParallel) { // if precice in parallel mode
      // code for serial write for parallel precice
      // maybe throw error saying inefficient?
    } else {
      // code for serial write of serial precice
    }
  }
}

void ExportVTKXML::processDataNamesAndDimensions
(
  mesh::Mesh& mesh)
{
  if (not _isDataNamesAndDimensionsProcessed) {
	  _isDataNamesAndDimensionsProcessed = true;
    _meshDimensions = mesh.getDimensions();
    if (_writeNormals) {
	  _vectorDataNames.push_back("VertexNormals ");
	  _vectorDataDimensions.push_back(3);
    }
    for (mesh::PtrData data : mesh.data()) {
	  int dataDimensions = data->getDimensions();
	  std::string dataName = data->getName();
	  if ( dataDimensions == 1) {
        _scalarDataNames.push_back(dataName);
	  } else if ( dataDimensions > 1 ) {
	    _vectorDataNames.push_back(dataName + " ");
	    _vectorDataDimensions.push_back(dataDimensions);
	  } else {
	    // TODO: Change to an assertion later
	    std::cout << "Error! Dimensions <= 0";
	    exit(-1);
	  }
    }
    if(mesh.edges().size() > 0 || mesh.triangles().size() > 0 || mesh.quads().size() > 0) {
    	_isCellPresent = true;
    }
  }
}

void ExportVTKXML::writeMasterFile
(
  const std::string& filename,
  mesh::Mesh&        mesh)
{
  std::ofstream outMasterFile;
  std::string fullMasterFilename(filename + "_master.pvtu");
  //utils::checkAppendExtension(fullMasterFilename, ".pvtu");

  outMasterFile.open(fullMasterFilename.c_str());
  preciceCheck(outMasterFile, "doExport()", "Could not open master file \"" << fullMasterFilename
      << "\" for VTKXML export!");

  outMasterFile << "<?xml version=\"1.0\"?>" << std::endl;
  outMasterFile << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"";
  outMasterFile << (utils::isMachineBigEndian() ? "BigEndian\">" : "LittleEndian\">")  << std::endl;
  outMasterFile << "   <PUnstructuredGrid GhostLevel=\"0\">" << std::endl;

  outMasterFile << "      <PPoints>" << std::endl;
  outMasterFile << "         <PDataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"" << _meshDimensions << "\"/>" << std::endl;
  outMasterFile << "      </PPoints>" << std::endl;

  if (_isCellPresent) {
    outMasterFile << "      <PCells>" << std::endl;
    outMasterFile << "         <PDataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\"/>" << std::endl;
    outMasterFile << "         <PDataArray type=\"Int32\" Name=\"offsets\"      NumberOfComponents=\"1\"/>" << std::endl;
    outMasterFile << "         <PDataArray type=\"UInt8\" Name=\"types\"        NumberOfComponents=\"1\"/>" << std::endl;
    outMasterFile << "      </PCells>" << std::endl;
  }

  // write scalar data names
  outMasterFile << "      <PPointData Scalars=\"";
  for (int i = 0; i < _scalarDataNames.size(); ++i) {
    outMasterFile << _scalarDataNames[i] << " ";
  }
  // write vector data names
  outMasterFile << "\" Vectors=\"";
  for (int i = 0; i < _vectorDataNames.size(); ++i) {
    outMasterFile << _vectorDataNames[i] << " ";
  }
  outMasterFile << "\">" << std::endl;

  for (int i = 0; i < _scalarDataNames.size(); ++i) {
    outMasterFile << "         <PDataArray type=\"Float32\" Name=\""<< _scalarDataNames[i] << "\" NumberOfComponents=\"" << 1 << "\"/>" << std::endl;
  }

  for (int i = 0; i < _vectorDataNames.size(); ++i) {
    outMasterFile << "         <PDataArray type=\"Float32\" Name=\""<< _vectorDataNames[i] << "\" NumberOfComponents=\"" << _vectorDataDimensions[i] << "\"/>" << std::endl;
  }
  outMasterFile << "      </PPointData>" << std::endl;

  for (int i = 0; i < utils::Parallel::getCommunicatorSize(); i++) {
    outMasterFile << "      <Piece Source=\"" << filename << "_r" << i << ".vtu\"/>" << std::endl;
  }

  outMasterFile << "   </PUnstructuredGrid>" << std::endl;
  outMasterFile << "</VTKFile>" << std::endl;

  outMasterFile.close();
}

void ExportVTKXML::writeSubFile
(
  const std::string& filename,
  mesh::Mesh&        mesh)
{
  int rank = utils::Parallel::getProcessRank(); // process rank
  int numPoints = mesh.vertices().size(); // number of vertices
  int numCells; // number of cells
  if (_meshDimensions == 2) {
    numCells = mesh.edges().size();
  } else {
	numCells = mesh.triangles().size();
  }

  std::ofstream outSubFile;
  std::string fullSubFilename(filename + "_r" + std::to_string(rank) + ".vtu");
  //utils::checkAppendExtension(fullMasterFilename, ".pvtu");

  outSubFile.open(fullSubFilename.c_str());
  preciceCheck(outSubFile, "doExport()", "Could not open master file \"" << fullSubFilename
       << "\" for VTKXML export!");

  outSubFile << "<?xml version=\"1.0\"?>" << std::endl;
  outSubFile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"";
  outSubFile << (utils::isMachineBigEndian() ? "BigEndian\">" : "LittleEndian\">")  << std::endl;

  outSubFile << "   <UnstructuredGrid>" << std::endl;
  outSubFile << "      <Piece NumberOfPoints=\"" << numPoints << "\" NumberOfCells=\"" << numCells << "\"> " << std::endl;
  outSubFile << "         <Points> " << std::endl;
  outSubFile << "            <DataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"" << _meshDimensions << "\" format=\"ascii\"> " << std::endl;
  for (mesh::Vertex& vertex : mesh.vertices()) {
    writeVertex(vertex.getCoords(), outSubFile);
  }
  outSubFile << "            </DataArray>" << std::endl;
  outSubFile << "         </Points> " << std::endl << std::endl;

  // Write geometry
  exportGeometry(outSubFile, mesh);

  // Write data
  exportData(outSubFile, mesh);

  outSubFile << "      </Piece>" << std::endl;
  outSubFile << "   </UnstructuredGrid> " << std::endl;
  outSubFile << "</VTKFile>" << std::endl;

  outSubFile.close();
}

void ExportVTKXML::exportGeometry
(
  std::ofstream& outFile,
  mesh::Mesh&    mesh)
{
  if (_isCellPresent) {
    if (_meshDimensions == 2) { // write edges as cells
      outFile << "         <Cells>" << std::endl;
	  outFile << "            <DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
	  for (mesh::Edge & edge : mesh.edges()) {
  	    writeLine(edge, outFile);
  	  }
	  outFile << std::endl;
	  outFile << "            </DataArray> " << std::endl;
	  outFile << "            <DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
	  for (int i = 1; i <= mesh.edges().size(); i++) {
		  outFile << 2*i << "  ";
	  }
	  outFile << std::endl;
	  outFile << "            </DataArray>" << std::endl;
	  outFile << "            <DataArray type=\"UInt8\"  Name=\"types\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
	  for (int i = 1; i <= mesh.edges().size(); i++) {
		  outFile << 3 << "  ";
	  }
	  outFile << std::endl;
	  outFile << "            </DataArray>" << std::endl;
	  outFile << "         </Cells>" << std::endl;

    } else { // write triangles and quads as cells

	  outFile << "         <Cells>" << std::endl;
	  outFile << "            <DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
	  for (mesh::Triangle& triangle : mesh.triangles()) {
	    writeTriangle(triangle, outFile);
	  }
	  for (mesh::Quad& quad : mesh.quads()) {
	    writeQuadrangle(quad, outFile);
	  }
	  outFile << std::endl;
	  outFile << "            </DataArray> " << std::endl;
	  outFile << "            <DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
	  for (int i = 1; i <= mesh.triangles().size(); i++) {
  	  outFile << 3*i << "  ";
	  }
	  for (int i = 1; i <= mesh.quads().size(); i++) {
	    outFile << 4*i << "  ";
	  }
	  outFile << std::endl;
	  outFile << "            </DataArray>" << std::endl;
	  outFile << "            <DataArray type=\"UInt8\"  Name=\"types\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
	  for (int i = 1; i <= mesh.triangles().size(); i++) {
        outFile << 5 << "  ";
	  }
	  for (int i = 1; i <= mesh.quads().size(); i++) {
	    outFile << 9 << "  ";
	  }
	  outFile << std::endl;
	  outFile << "            </DataArray>" << std::endl;
	  outFile << "         </Cells>" << std::endl;
    }
  }
}

void ExportVTKXML:: exportData
(
  std::ofstream& outFile,
  mesh::Mesh&    mesh)
{
  outFile << "         <PointData Scalars=\"";
  for (int i = 0; i < _scalarDataNames.size(); i++) {
	outFile << _scalarDataNames[i] << " ";
  }
  outFile << "\" Vectors=\"";
  for (int i = 0; i < _vectorDataNames.size(); i++) {
	outFile << _vectorDataNames[i] << " ";
  }
  outFile << "\">" << std::endl;

  for (mesh::PtrData data : mesh.data()) { // Plot vertex data
	Eigen::VectorXd& values = data->values();
	int dataDimensions = data->getDimensions();
	std::string dataName(data->getName());
	outFile << "            <DataArray type=\"Float32\" Name=\"" << dataName << "\" NumberOfComponents=\"" << dataDimensions << "\" format=\"ascii\">" << std::endl;
	outFile << "               ";
	if(dataDimensions > 1) {
	  utils::DynVector viewTemp(dataDimensions);
	  for (mesh::Vertex& vertex : mesh.vertices()) {
		int offset = vertex.getID() * dataDimensions;
		for(int i=0; i < dataDimensions; i++){
		  viewTemp[i] = values(offset + i);
		}
		for(int i = 0; i < dataDimensions; i++){
	      outFile << viewTemp[i] << " ";
		}
		outFile << " ";
	  }
	} else if(dataDimensions == 1) {
	  for (mesh::Vertex& vertex : mesh.vertices()) {
		outFile << values(vertex.getID()) << " ";
	  }
	}
	outFile << std::endl << "            </DataArray>" << std::endl;
  }

  outFile << "         </PointData> " << std::endl;
}

void ExportVTKXML::writeVertex
(
  const utils::DynVector& position,
  std::ofstream&           outFile)
{
  outFile << "               ";
  for (int i = 0; i < position.size(); i++){
	outFile << position(i) << "  ";
  }
  outFile << std::endl;
}


void ExportVTKXML:: writeTriangle
(
  mesh::Triangle& triangle,
  std::ofstream&  outFile)
{
  outFile << triangle.vertex(0).getID() << "  ";
  outFile << triangle.vertex(1).getID() << "  ";
  outFile << triangle.vertex(2).getID() << "  ";

}

void ExportVTKXML:: writeQuadrangle
(
  mesh::Quad&    quad,
  std::ofstream& outFile)
{
  outFile << quad.vertex(0).getID() << "  ";
  outFile << quad.vertex(1).getID() << "  ";
  outFile << quad.vertex(2).getID() << "  ";
  outFile << quad.vertex(3).getID() << "  ";
}

void ExportVTKXML:: writeLine
(
  mesh::Edge& edge,
  std::ofstream& outFile)
{
  outFile << edge.vertex(0).getID() << "  ";
  outFile << edge.vertex(1).getID() << "  ";
}

}} // namespace precice, io

