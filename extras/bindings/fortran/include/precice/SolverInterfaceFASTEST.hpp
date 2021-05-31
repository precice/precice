#pragma once

/**
 * @file
 * This file contains a specific FASTEST Fortran 77 compatible interface
 * written in C/C++. The specialty is that we offer here 2 interfaces. We need
 * this since FASTEST uses an internal subcycling, coupling the FASTEST fluid
 * solver to a structure solver at big timesteps and the FASTEST acoustic solver
 * to an external pure acoustic solver at small timesteps. So, we need these
 * bindings to handle two participants, FASTEST-Acoustic and FASTEST-Fluid,
 * creating two C++ SolverInterfaces.
 *
 * As for the "normal" Fortran bindings, every method has a Fortran syntax
 * equivalent in the method comment, and a listing for input and output
 * variables. A variable can be input and output at the same time.
 */

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief See precice::SolverInterface::SolverInterface() and ...::configure().
 *
 * Fortran syntax:
 * precicef_create(
 *   CHARACTER participantNameAcoustic(*),
 *   INTEGER   isAcousticUsed,
 *   CHARACTER participantNameFluid(*),
 *   INTEGER   isFluidUsed,
 *   CHARACTER configFileName(*),
 *   INTEGER   solverProcessIndex,
 *   INTEGER   solverProcessSize )
 *
 * IN:  participantNameA
 * IN:  isAcousticUsed: Acoustic solver is used (0 or 1)
 * IN:  participantNameF
 * IN:  isFluidUsed: Fluid solver is used (0 or 1)
 * IN:  configFileName
 * IN:  solverProcessIndex
 * IN:  solverProcessSize
 * OUT: -
 */
void precice_fastest_create_(
    const char *participantNameAcoustic,
    const int * isAcousticUsed,
    const char *participantNameFluid,
    const int * isFluidUsed,
    const char *configFileName,
    const int * solverProcessIndex,
    const int * solverProcessSize,
    int         lengthAccessorNameAcoustic,
    int         lengthAccessorNameFluid,
    int         lengthConfigFileName);

/**
 * @brief See precice::SolverInterface::initialize().
 *
 * Fortran syntax:
 * precicef_initialize( DOUBLE PRECISION timstepLengthLimit, INTEGER useFluid )
 *
 * IN: useFluid (1: use the Fluid interface , 0: use the Acoustic interface)
 * OUT: timestepLengthLimit
 */
void precice_fastest_initialize_(double *timestepLengthLimit, const int *useFluid);

/**
 * @brief See precice::SolverInterface::initializeData().
 *
 * Fortran syntax:
 * precicef_intialize_data(INTEGER useFluid)
 *
 * IN: useFluid (1: use the Fluid interface , 0: use the Acoustic interface)
 * OUT: -
 */
void precice_fastest_initialize_data_(const int *useFluid);

/**
 * @brief See precice::SolverInterface::advance().
 *
 * Fortran syntax:
 * precicef_advance( DOUBLE PRECISION timstepLengthLimit, INTEGER useFluid )
 *
 * IN:  timestepLengthLimit
 * IN:  useFluid (1: use the Fluid interface , 0: use the Acoustic interface)
 * OUT: timestepLengthLimit
 */
void precice_fastest_advance_(double *timestepLengthLimit, const int *useFluid);

/**
 * @brief See precice::SolverInterface::finalize().
 *
 * Fortran syntax:
 * precicef_finalize(INTEGER useFluid);
 *
 * IN:  useFluid (1: use the Fluid interface , 0: use the Acoustic interface)
 */
void precice_fastest_finalize_(const int *useFluid);

/**
 * @brief See precice::SolverInterface::isActionRequired().
 *
 * Fortran syntax:
 * precicef_action_required(
 *   CHARACTER action(*),
 *   INTEGER   isRequired,
 *   INTEGER   useFluid )
 *
 * IN:  action
 * IN:  useFluid (1: use the Fluid interface , 0: use the Acoustic interface)
 * OUT: isRequired(1:true, 0:false)
 */
void precice_fastest_action_required_(
    const char *action,
    int *       isRequired,
    const int * useFluid,
    int         lengthAction);

/**
 * @brief See precice::SolverInterface::markActionFulfilled().
 *
 * Fortran syntax:
 * precicef_mark_action_fulfilled_( CHARACTER action(*), INTEGER useFluid )
 *
 * IN:  action
 * IN:  useFluid (1: use the Fluid interface , 0: use the Acoustic interface)
 * OUT: -
 */
void precice_fastest_mark_action_fulfilled_(
    const char *action,
    const int * useFluid,
    int         lengthAction);

/**
 * @brief See precice::SolverInterface::getMeshID().
 *
 * Fortran syntax:
 * precicef_get_mesh_id(
 *   CHARACTER meshName(*),
 *   INTEGER   meshID,
 *   INTEGER   useFluid )
 *
 * IN:  meshName
 * IN:  useFluid (1: use the Fluid interface , 0: use the Acoustic interface)
 * OUT: meshID
 */
void precice_fastest_get_mesh_id_(
    const char *meshName,
    int *       meshID,
    const int * useFluid,
    int         lengthMeshName);

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
 *   INTEGER   useFluid)
 *
 * IN:  dataName
 * IN:  meshID
 * IN:  useFluid (1: use the Fluid interface , 0: use the Acoustic interface)
 * OUT: dataID
 */
void precice_fastest_get_data_id_(
    const char *dataName,
    const int * meshID,
    int *       dataID,
    const int * useFluid,
    int         lengthDataName);

/**
 * @brief See precice::SolverInterface::setMeshVertex().
 *
 * Fortran syntax:
 * precicef_set_vertex(
 *   INTEGER          meshID,
 *   DOUBLE PRECISION position(dim),
 *   INTEGER          vertexID,
 *   INTEGER          useFluid )
 *
 * IN:  meshID, position
 * IN:  useFluid (1: use the Fluid interface , 0: use the Acoustic interface)
 * OUT: vertexID
 */
void precice_fastest_set_vertex_(
    const int *   meshID,
    const double *position,
    int *         vertexID,
    const int *   useFluid);

/**
 * @brief See precice::SolverInterface::setMeshVertices().
 *
 * Fortran syntax:
 * precicef_set_vertices(
 *   INTEGER          meshID,
 *   INTEGER          size,
 *   DOUBLE PRECISION positions(dim*size),
 *   INTEGER          positionIDs(size),
 *   INTEGER          useFluid )
 *
 * IN:  meshID, size, positions
 * IN:  useFluid (1: use the Fluid interface , 0: use the Acoustic interface)
 * OUT: positionIDs
 */
void precice_fastest_set_vertices_(
    const int *meshID,
    const int *size,
    double *   positions,
    int *      positionIDs,
    const int *useFluid);

/**
 * @brief See precice::SolverInterface::writeBlockVectorData.
 *
 * Fortran syntax:
 * precicef_write_bvdata(
 *   INTEGER dataID,
 *   INTEGER size,
 *   INTEGER valueIndices,
 *   DOUBLE PRECISION values(dim*size),
 *   INTEGER useFluid )
 *
 * IN:  dataID, size, valueIndices, values
 * IN:  useFluid (1: use the Fluid interface , 0: use the Acoustic interface)
 * OUT: -
 */
void precice_fastest_write_bvdata_(
    const int *dataID,
    const int *size,
    int *      valueIndices,
    double *   values,
    const int *useFluid);

/**
 * @brief precice::SolverInterface::writeVectorData.
 *
 * Fortran syntax:
 * precicef_write_vdata(
 *   INTEGER dataID,
 *   INTEGER valueIndex,
 *   DOUBLE PRECISION dataValue(dim),
 *   INTEGER useFluid )
 *
 * IN:  dataID, valueIndex, dataValue
 * IN:  useFluid (1: use the Fluid interface , 0: use the Acoustic interface)
 * OUT: -
 */
void precice_fastest_write_vdata_(
    const int *   dataID,
    const int *   valueIndex,
    const double *dataValue,
    const int *   useFluid);

/**
 * @brief See precice::SolverInterface::writeBlockScalarData.
 *
 * Fortran syntax:
 * precicef_write_bsdata(
 *   INTEGER dataID,
 *   INTEGER size,
 *   INTEGER valueIndices,
 *   DOUBLE PRECISION values(size),
 *   INTEGER useFluid )
 *
 * IN:  dataID, size, valueIndices, values
 * IN:  useFluid (1: use the Fluid interface , 0: use the Acoustic interface)
 * OUT: -
 */
void precice_fastest_write_bsdata_(
    const int *dataID,
    const int *size,
    int *      valueIndices,
    double *   values,
    const int *useFluid);

/**
 * @brief precice::SolverInterface::writeScalarData.
 *
 * Fortran syntax:
 * precicef_write_sdata(
 *   INTEGER dataID,
 *   INTEGER valueIndex,
 *   DOUBLE PRECISION dataValue,
 *   INTEGER useFluid )
 *
 * IN:  dataID, valueIndex, dataValue
 * IN:  useFluid (1: use the Fluid interface , 0: use the Acoustic interface)
 * OUT: -
 */
void precice_fastest_write_sdata_(
    const int *   dataID,
    const int *   valueIndex,
    const double *dataValue,
    const int *   useFluid);

/**
 * @brief See precice::SolverInterface::readBlockVectorData.
 *
 * Fortran syntax:
 * precicef_read_bvdata(
 *   INTEGER dataID,
 *   INTEGER size,
 *   INTEGER valueIndices,
 *   DOUBLE PRECISION values(dim*size),
 *   INTEGER useFluid )
 *
 * IN:  dataID, size, valueIndices
 * IN:  useFluid (1: use the Fluid interface , 0: use the Acoustic interface)
 * OUT: values
 */
void precice_fastest_read_bvdata_(
    const int *dataID,
    const int *size,
    int *      valueIndices,
    double *   values,
    const int *useFluid);

/**
 * @brief precice::SolverInterface::readVectorData.
 *
 * Fortran syntax:
 * precicef_read_vdata(
 *   INTEGER dataID,
 *   INTEGER valueIndex,
 *   DOUBLE PRECISION dataValue(dim),
 *   INTEGER useFluid )
 *
 * IN:  dataID, valueIndex
 * IN:  useFluid (1: use the Fluid interface , 0: use the Acoustic interface)
 * OUT: dataValue
 */
void precice_fastest_read_vdata_(
    const int *dataID,
    const int *valueIndex,
    double *   dataValue,
    const int *useFluid);

/**
 * @brief See precice::SolverInterface::readBlockScalarData.
 *
 * Fortran syntax:
 * precicef_read_bsdata(
 *   INTEGER dataID,
 *   INTEGER size,
 *   INTEGER valueIndices,
 *   DOUBLE PRECISION values(size),
 *   INTEGER useFluid )
 *
 * IN:  dataID, size, valueIndices
 * IN:  useFluid (1: use the Fluid interface , 0: use the Acoustic interface)
 * OUT: values
 */
void precice_fastest_read_bsdata_(
    const int *dataID,
    const int *size,
    int *      valueIndices,
    double *   values,
    const int *useFluid);

/**
 * @brief precice::SolverInterface::readScalarData.
 *
 * Fortran syntax:
 * precicef_read_sdata(
 *   INTEGER dataID,
 *   INTEGER valueIndex,
 *   DOUBLE PRECISION dataValue,
 *   INTEGER useFluid )
 *
 * IN:  dataID, valueIndex
 * IN:  useFluid (1: use the Fluid interface , 0: use the Acoustic interface)
 * OUT: dataValue
 */
void precice_fastest_read_sdata_(
    const int *dataID,
    const int *valueIndex,
    double *   dataValue,
    const int *useFluid);

#ifdef __cplusplus
}
#endif
