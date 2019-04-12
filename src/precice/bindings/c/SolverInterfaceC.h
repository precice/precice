#pragma once

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Creates the coupling interface and confiures it.
 *
 * Has to be called before any other method of this interface.
 *
 * @param[in] participantName Name of the participant accessing the interface. Has to
 *                          match one of the names specified in the
 *                          configuration xml file.
 * @param[in] configFileName (Path and) name of the xml configuration file
 *                            containing the precice configuration.
 */
void precicec_createSolverInterface (
  const char* participantName,
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
 * @brief Initializes coupling data.
 */
void precicec_initialize_data();

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
 * @brief Returns id belonging to the given mesh name
 */
int precicec_getMeshID ( const char* meshName );

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
  int        meshID,
  int        size,
  const int* ids,
  double*    positions );

void precicec_setMeshVertices (
  int           meshID,
  int           size,
  const double* positions,
  int*          ids );

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
 * @param[in] dataID ID of the data to be written.
 * @param[in] size Number of indices, and number of values * dimensions.
 * @param[in] values Values of the data to be written.
 */
void precicec_writeBlockVectorData (
  int           dataID,
  int           size,
  const int*    valueIndices,
  const double* values );

/**
 * @brief Writes vectorial foating point data to the coupling mesh.
 *
 * @param[in] dataID ID of the data to be written. Obtained by getDataID().
 * @param[in] dataPosition Spatial position of the data to be written.
 * @param[in] dataValue Vectorial data value to be written.
 */
void precicec_writeVectorData (
  int           dataID,
  int           valueIndex,
  const double* dataValue );

/**
 * @brief See precice::SolverInterface::writeBlockScalarData().
 */
void precicec_writeBlockScalarData (
  int           dataID,
  int           size,
  const int*    valueIndices,
  const double* values );

/**
 * @brief Writes scalar floating point data to the coupling mesh.
 *
 * @param[in] dataID ID of the data to be written. Obtained by getDataID().
 * @param[in] dataPosition Spatial position of the data to be written.
 * @param[in] dataValue Scalar data value to be written.
 */
void precicec_writeScalarData (
  int           dataID,
  int           valueIndex,
  const double& dataValue );

/**
 * @brief Reads vector data values given as block.
 *
 * The block contains the vector values in the following form:
 * values = (d0x, d0y, d0z, d1x, d1y, d1z, ...., dnx, dny, dnz), where n is
 * the number of vector values. In 2D, the z-components are removed.
 *
 * @param[in] dataID ID of the data to be read.
 * @param[in] size  Number of indices, and number of values * dimensions.
 * @param[in] valueIndices Indices (from setReadPosition()) of data values.
 * @param[in] values Values of the data to be read.
 */
void precicec_readBlockVectorData (
  int        dataID,
  int        size,
  const int* valueIndices,
  double*    values );

/**
 * @brief Reads vectorial foating point data from the coupling mesh.
 *
 * @param[in] dataID ID of the data to be read. Obtained by getDataID().
 * @param[in] dataPosition Position where the read data should be mapped to.
 * @param[out] dataValue Vectorial data value read.
 */
void precicec_readVectorData (
  int     dataID,
  int     valueIndex,
  double* dataValue );

/**
 * @brief See precice::SolverInterface::readBlockScalarData().
 */
void precicec_readBlockScalarData (
  int        dataID,
  int        size,
  const int* valueIndices,
  double*    values );

/**
 * @brief Reads scalar foating point data from the coupling mesh.
 *
 * @param[in] dataID ID of the data to be read. Obtained by getDataID().
 * @param[in] dataPosition Position where the read data should be mapped to.
 * @param[out] dataValue Scalar data value read.
 */
void precicec_readScalarData (
  int     dataID,
  int     valueIndex,
  double* dataValue );

/**
 * @brief Computes and maps all write data mapped from mesh with given ID.
 */
void precicec_mapWriteDataFrom ( int fromMeshID );

/**
 * @brief Computes and maps all read data mapped to mesh with given ID.
 */
void precicec_mapReadDataTo ( int toMeshID );


#ifdef __cplusplus
}
#endif
