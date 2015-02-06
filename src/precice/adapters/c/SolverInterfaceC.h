/* Copyright (C) 2011 Technische Universitaet Muenchen
 * This file is part of the preCICE project. For conditions of distribution and
 * use, please see the license notice at http://www5.in.tum.de/wiki/index.php/precice_c_License */
#ifndef PRECICE_ADAPTERS_C_SOLVERINTERFACEC_H_
#define PRECICE_ADAPTERS_C_SOLVERINTERFACEC_H_

/**
 * @brief Creates the coupling interface and confiures it.
 *
 * Has to be called before any other method of this interface.
 *
 * @param accessorName [IN] Name of the solver accessing the interface. Has to
 *                          match one of the names specified in the
 *                          configuration xml file.
 * @param configFileName [IN] (Path and) name of the xml configuration file
 *                            containing the precice configuration.
 */
void precicec_createSolverInterface (
  const char* accessorName,
  const char* configFileName,
  int         solverProcessIndex,
  int         solverProcessSize );

/**
 * @brief Initiates the coupling to the coupling supervisor.
 *
 * @return Maximal length of first timestep to be computed by solver.
 */
double precicec_initialize();

/**
 * @brief Exchanges data between solver and coupling supervisor.
 *
 * @param computedTimestepLength [IN] Length of timestep computed by solver.
 * @return Maximal length of next timestep to be computed by solver.
 */
double precicec_advance ( double computedTimestepLength );

/**
 * @brief Finalizes the coupling to the coupling supervisor.
 */
void precicec_finalize();

/**
 * @brief Returns the number of spatial configurations for the coupling.
 */
int precicec_getDimensions();

/**
 * @brief Returns true (->1), if the coupled simulation is ongoing
 */
int precicec_isCouplingOngoing();

/**
 * @brief Returns true (->1), if the coupling timestep is completed.
 */
int precicec_isCouplingTimestepComplete();

int precicec_isWriteDataRequired ( double computedTimestepLength );

/**
 * @brief Returns true (->1), if new data to read is available.
 */
int precicec_isReadDataAvailable();


int precicec_isActionRequired ( const char* action );

void precicec_fulfilledAction ( const char* action );

/**
 * @brief Returns geometry id belonging to the given geometry name
 */
int precicec_getMeshID ( const char* geometryName );

/**
 * @brief Returns true (!=0), if data with given name is available.
 */
int precicec_hasData ( const char* dataName, int meshID );

/**
 * @brief Returns the data id belonging to the given name.
 *
 * The given name (dataName) has to be one of the names specified in the
 * configuration file. The data id obtained can be used to read and write
 * data to and from the coupling mesh.
 */
int precicec_getDataID ( const char* dataName, int meshID );

int precicec_setMeshVertex (
  int           meshID,
  const double* position );

void precicec_getMeshVertices (
  int     meshID,
  int     size,
  int*    ids,
  double* positions );

int precicec_getMeshVertexSize ( int meshID );


int precicec_setMeshEdge (
  int meshID,
  int firstVertexID,
  int secondVertexID );

void precicec_setMeshTriangle (
  int meshID,
  int firstEdgeID,
  int secondEdgeID,
  int thirdEdgeID );

/**
 * @brief Sets a triangle from vertex IDs. Creates missing edges.
 */
void precicec_setMeshTriangleWithEdges (
  int meshID,
  int firstVertexID,
  int secondVertexID,
  int thirdVertexID );

/**
 * @brief Writes vector data values given as block.
 *
 * The block must contain the vector values in the following form:
 * values = (d0x, d0y, d0z, d1x, d1y, d1z, ...., dnx, dny, dnz), where n is
 * the number of vector values. In 2D, the z-components are removed.
 *
 * @param dataID [IN] ID of the data to be written.
 * @param size [IN] Number of indices, and number of values * dimensions.
 * @param values [IN] Values of the data to be written.
 */
void precicec_writeBlockVectorData (
  int     dataID,
  int     size,
  int*    valueIndices,
  double* values );

/**
 * @brief Writes vectorial foating point data to the coupling mesh.
 *
 * @param dataID [IN] ID of the data to be written. Obtained by getDataID().
 * @param dataPosition [IN] Spatial position of the data to be written.
 * @param dataValue [IN] Vectorial data value to be written.
 */
void precicec_writeVectorData (
  int           dataID,
  int           valueIndex,
  const double* dataValue );

/**
 * @brief See precice::SolverInterface::writeBlockScalarData().
 */
void precicec_writeBlockScalarData (
  int     dataID,
  int     size,
  int*    valueIndices,
  double* values );

/**
 * @brief Writes scalar floating point data to the coupling mesh.
 *
 * @param dataID [IN] ID of the data to be written. Obtained by getDataID().
 * @param dataPosition [IN] Spatial position of the data to be written.
 * @param dataValue [IN] Scalar data value to be written.
 */
void precicec_writeScalarData (
  int    dataID,
  int    valueIndex,
  double dataValue );

/**
 * @brief Reads vector data values given as block.
 *
 * The block contains the vector values in the following form:
 * values = (d0x, d0y, d0z, d1x, d1y, d1z, ...., dnx, dny, dnz), where n is
 * the number of vector values. In 2D, the z-components are removed.
 *
 * @param dataID [IN] ID of the data to be read.
 * @param size [IN] Number of indices, and number of values * dimensions.
 * @param valueIndices [IN] Indices (from setReadPosition()) of data values.
 * @param values [IN] Values of the data to be read.
 */
void precicec_readBlockVectorData (
  int     dataID,
  int     size,
  int*    valueIndices,
  double* values );

/**
 * @brief Reads vectorial foating point data from the coupling mesh.
 *
 * @param dataID [IN] ID of the data to be read. Obtained by getDataID().
 * @param dataPosition [IN] Position where the read data should be mapped to.
 * @param dataValue [OUT] Vectorial data value read.
 */
void precicec_readVectorData (
  int     dataID,
  int     valueIndex,
  double* dataValue );

/**
 * @brief See precice::SolverInterface::readBlockScalarData().
 */
void precicec_readBlockScalarData (
  int     dataID,
  int     size,
  int*    valueIndices,
  double* values );

/**
 * @brief Reads scalar foating point data from the coupling mesh.
 *
 * @param dataID [IN] ID of the data to be read. Obtained by getDataID().
 * @param dataPosition [IN] Position where the read data should be mapped to.
 * @param dataValue [OUT] Scalar data value read.
 */
void precicec_readScalarData (
  int     dataID,
  int     valueIndex,
  double* dataValue );

/**
 * @brief Computes and maps all write data mapped from mesh with given ID.
 *
 */
void precicec_mapWriteDataFrom ( int fromMeshID );

/**
 * @brief Computes and maps all read data mapped to mesh with given ID.
 */
void precicec_mapReadDataTo ( int toMeshID );

/**
 * @brief Exports the coupling mesh to a vtk file
 *
 * @param filenameSuffix [IN] Suffix added to the vtk filename.
 */
void precicec_exportMesh ( const char* filenameSuffix );

#endif /* PRECICE_ADAPTERS_C_SOLVERINTERFACEC_H_ */
