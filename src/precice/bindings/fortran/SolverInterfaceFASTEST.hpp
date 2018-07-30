#pragma once

/**
 * @file This file contains a specific FASTEST Fortran 77 compatible interface written in C/C++.
 * The specialty is that we offer here 2 interfaces. We needed this since FASTEST uses an internal subcycling,
 * coupling to a structure solver at big timesteps and to a pure acoustic solver at small timesteps
 *
 * It has been tested with: gfortran 4.4.3, ifort 12.1.0
 *
 * Every method has a Fortran syntax equivalent in the method comment, and a
 * listing for input and output variables. A variable can be input and output
 * at the same time.
 */

#ifdef __cplusplus
extern"C" {
#endif

/**
 * @brief See precice::SolverInterface::SolverInterface() and ...::configure().
 *
 * Fortran syntax:
 * precicef_create(
 *   CHARACTER participantNameA(*),
 *   CHARACTER participantNameF(*),
 *   CHARACTER configFileName(*),
 *   INTEGER   solverProcessIndex,
 *   INTEGER   solverProcessSize )
 *
 * IN:  participantNameA, participantNameF, configFileName, solverProcessIndex, solverProcessSize
 * OUT: -
 */
void precice_fastest_create_(
  const char* participantNameA,
  const char* participantNameF,
  const char* configFileName,
  const int*  solverProcessIndex,
  const int*  solverProcessSize,
  int   lengthAccessorNameA,
  int   lengthAccessorNameF,
  int   lengthConfigFileName );

/**
 * @brief See precice::SolverInterface::initialize().
 *
 * Fortran syntax:
 * precicef_initialize( DOUBLE PRECISION timstepLengthLimit, INTEGER useF )
 *
 * IN: useF (1: this is the F interface , 0: this is the A interface)
 * OUT: timestepLengthLimit
 */
void precice_fastest_initialize_( double* timestepLengthLimit, const int*  useF );

/**
 * @brief See precice::SolverInterface::initializeData().
 *
 * Fortran syntax:
 * precicef_intialize_data(INTEGER useF)
 *
 * IN: useF (1: this is the F interface , 0: this is the A interface)
 * OUT: -
 */
void precice_fastest_initialize_data_(const int*  useF);

/**
 * @brief See precice::SolverInterface::advance().
 *
 * Fortran syntax:
 * precicef_advance( DOUBLE PRECISION timstepLengthLimit, INTEGER useF )
 *
 * IN:  timestepLengthLimit
 * IN:  useF (1: this is the F interface , 0: this is the A interface)
 * OUT: timestepLengthLimit
 */
void precice_fastest_advance_( double* timestepLengthLimit, const int*  useF );


/**
 * @brief See precice::SolverInterface::finalize().
 *
 * Fortran syntax:
 * precicef_finalize(INTEGER useF);
 *
 * IN:  useF (1: this is the F interface , 0: this is the A interface)
 */
void precice_fastest_finalize_(const int*  useF);

/**
 * @brief See precice::SolverInterface::isActionRequired().
 *
 * Fortran syntax:
 * precicef_action_required(
 *   CHARACTER action(*),
 *   INTEGER   isRequired,
 *   INTEGER   useF )
 *
 * IN:  action
 * IN:  useF (1: this is the F interface , 0: this is the A interface)
 * OUT: isRequired(1:true, 0:false)
 */
void precice_fastest_action_required_(
  const char* action,
  int*        isRequired,
  const int*  useF,
  int         lengthAction );

/**
 * @brief See precice::SolverInterface::fulfilledAction().
 *
 * Fortran syntax:
 * precicef_fulfilled_action( CHARACTER action(*), INTEGER useF )
 *
 * IN:  action
 * IN:  useF (1: this is the F interface , 0: this is the A interface)
 * OUT: -
 */
void precice_fastest_fulfilled_action_(
  const char* action,
  const int*  useF,
  int         lengthAction );

/**
 * @brief See precice::SolverInterface::getMeshID().
 *
 * Fortran syntax:
 * precicef_get_mesh_id(
 *   CHARACTER meshName(*),
 *   INTEGER   meshID,
 *   INTEGER   useF )
 *
 * IN:  meshName
 * IN:  useF (1: this is the F interface , 0: this is the A interface)
 * OUT: meshID
 */
void precice_fastest_get_mesh_id_(
  const char* meshName,
  int*        meshID,
  const int*  useF,
  int         lengthMeshName );


/**
 * @brief See precice::SolverInterface::getDataID().
 *
 * The given name (dataName) has to be one of the names specified in the
 * configuration file. The data id obtained can be used to read and write
 * data to and from the coupling mesh.
 *
 * Fortran syntax:
 * precicef_get_data_id(
 *   CHARACTER dataName(*),
 *   INTEGER   meshID,
 *   INTEGER   dataID,
 *   INTEGER   useF)
 *
 * IN:  dataName
 * IN:  meshID
 * IN:  useF (1: this is the F interface , 0: this is the A interface)
 * OUT: dataID
 */
void precice_fastest_get_data_id_(
  const char* dataName,
  const int*  meshID,
  int*        dataID,
  const int*  useF,
  int         lengthDataName);

/**
 * @brief See precice::SolverInterface::setMeshVertex().
 *
 * Fortran syntax:
 * precicef_set_vertex(
 *   INTEGER          meshID,
 *   DOUBLE PRECISION position(dim),
 *   INTEGER          vertexID,
 *   INTEGER          useF )
 *
 * IN:  meshID, position
 * IN:  useF (1: this is the F interface , 0: this is the A interface)
 * OUT: vertexID
 */
void precice_fastest_set_vertex_(
  const int*    meshID,
  const double* position,
  int*          vertexID,
  const int*    useF);

/**
 * @brief See precice::SolverInterface::setMeshVertices().
 *
 * Fortran syntax:
 * precicef_set_vertices(
 *   INTEGER          meshID,
 *   INTEGER          size,
 *   DOUBLE PRECISION positions(dim*size),
 *   INTEGER          positionIDs(size),
 *   INTEGER          useF )
 *
 * IN:  meshID, size, positions
 * IN:  useF (1: this is the F interface , 0: this is the A interface)
 * OUT: positionIDs
 */
void precice_fastest_set_vertices_(
  const int*    meshID,
  const int*    size,
  double*       positions,
  int*          positionIDs,
  const int*    useF);


/**
 * @brief See precice::SolverInterface::writeBlockVectorData.
 *
 * Fortran syntax:
 * precicef_write_bvdata(
 *   INTEGER dataID,
 *   INTEGER size,
 *   INTEGER valueIndices,
 *   DOUBLE PRECISION values(dim*size),
 *   INTEGER useF )
 *
 * IN:  dataID, size, valueIndices, values
 * IN:  useF (1: this is the F interface , 0: this is the A interface)
 * OUT: -
 */
void precice_fastest_write_bvdata_(
  const int* dataID,
  const int* size,
  int*       valueIndices,
  double*    values,
  const int* useF);

/**
 * @brief precice::SolverInterface::writeVectorData.
 *
 * Fortran syntax:
 * precicef_write_vdata(
 *   INTEGER dataID,
 *   INTEGER valueIndex,
 *   DOUBLE PRECISION dataValue(dim),
 *   INTEGER useF )
 *
 * IN:  dataID, valueIndex, dataValue
 * IN:  useF (1: this is the F interface , 0: this is the A interface)
 * OUT: -
 */
void precice_fastest_write_vdata_(
  const int*    dataID,
  const int*    valueIndex,
  const double* dataValue,
  const int*    useF);

/**
 * @brief See precice::SolverInterface::writeBlockScalarData.
 *
 * Fortran syntax:
 * precicef_write_bsdata(
 *   INTEGER dataID,
 *   INTEGER size,
 *   INTEGER valueIndices,
 *   DOUBLE PRECISION values(size),
 *   INTEGER useF )
 *
 * IN:  dataID, size, valueIndices, values
 * IN:  useF (1: this is the F interface , 0: this is the A interface)
 * OUT: -
 */
void precice_fastest_write_bsdata_(
  const int* dataID,
  const int* size,
  int*       valueIndices,
  double*    values,
  const int* useF);

/**
 * @brief precice::SolverInterface::writeScalarData.
 *
 * Fortran syntax:
 * precicef_write_sdata(
 *   INTEGER dataID,
 *   INTEGER valueIndex,
 *   DOUBLE PRECISION dataValue,
 *   INTEGER useF )
 *
 * IN:  dataID, valueIndex, dataValue
 * IN:  useF (1: this is the F interface , 0: this is the A interface)
 * OUT: -
 */
void precice_fastest_write_sdata_(
  const int*    dataID,
  const int*    valueIndex,
  const double* dataValue,
  const int*    useF);

/**
 * @brief See precice::SolverInterface::readBlockVectorData.
 *
 * Fortran syntax:
 * precicef_read_bvdata(
 *   INTEGER dataID,
 *   INTEGER size,
 *   INTEGER valueIndices,
 *   DOUBLE PRECISION values(dim*size),
 *   INTEGER useF )
 *
 * IN:  dataID, size, valueIndices
 * IN:  useF (1: this is the F interface , 0: this is the A interface)
 * OUT: values
 */
void precice_fastest_read_bvdata_(
  const int* dataID,
  const int* size,
  int*       valueIndices,
  double*    values,
  const int* useF);

/**
 * @brief precice::SolverInterface::readVectorData.
 *
 * Fortran syntax:
 * precicef_read_vdata(
 *   INTEGER dataID,
 *   INTEGER valueIndex,
 *   DOUBLE PRECISION dataValue(dim),
 *   INTEGER useF )
 *
 * IN:  dataID, valueIndex
 * IN:  useF (1: this is the F interface , 0: this is the A interface)
 * OUT: dataValue
 */
void precice_fastest_read_vdata_(
  const int* dataID,
  const int* valueIndex,
  double*    dataValue,
  const int* useF);

/**
 * @brief See precice::SolverInterface::readBlockScalarData.
 *
 * Fortran syntax:
 * precicef_read_bsdata(
 *   INTEGER dataID,
 *   INTEGER size,
 *   INTEGER valueIndices,
 *   DOUBLE PRECISION values(size),
 *   INTEGER useF )
 *
 * IN:  dataID, size, valueIndices
 * IN:  useF (1: this is the F interface , 0: this is the A interface)
 * OUT: values
 */
void precice_fastest_read_bsdata_(
  const int* dataID,
  const int* size,
  int*       valueIndices,
  double*    values,
  const int* useF);

/**
 * @brief precice::SolverInterface::readScalarData.
 *
 * Fortran syntax:
 * precicef_read_sdata(
 *   INTEGER dataID,
 *   INTEGER valueIndex,
 *   DOUBLE PRECISION dataValue,
 *   INTEGER useF )
 *
 * IN:  dataID, valueIndex
 * IN:  useF (1: this is the F interface , 0: this is the A interface)
 * OUT: dataValue
 */
void precice_fastest_read_sdata_(
  const int* dataID,
  const int* valueIndex,
  double*    dataValue,
  const int* useF);


#ifdef __cplusplus
}
#endif
