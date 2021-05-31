#include "precice/SolverInterfaceFASTEST.hpp"
#include <iostream>
#include <memory>
#include <string>
#include "logging/LogMacros.hpp"
#include "logging/Logger.hpp"
#include "precice/SolverInterface.hpp"

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
#endif

using namespace std;

static std::unique_ptr<precice::SolverInterface> implAcoustic = nullptr;
static std::unique_ptr<precice::SolverInterface> implFluid    = nullptr;

static precice::logging::Logger _log("SolverInterfaceFASTEST");

namespace precice {
namespace impl {
/**
     * @brief Returns length of string without trailing whitespaces.
     */
int  strippedLength(const char *string, int length);
void checkCorrectUsage(int useFluid);
} // namespace impl
} // namespace precice

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
    int         lengthConfigFileName)
{
  int    strippedLength = precice::impl::strippedLength(participantNameAcoustic, lengthAccessorNameAcoustic);
  string stringAccessorNameAcoustic(participantNameAcoustic, strippedLength);
  strippedLength = precice::impl::strippedLength(participantNameFluid, lengthAccessorNameFluid);
  string stringAccessorNameFluid(participantNameFluid, strippedLength);
  strippedLength = precice::impl::strippedLength(configFileName, lengthConfigFileName);
  string stringConfigFileName(configFileName, strippedLength);
  //cout << "Accessor: " << stringAccessorName << "!" << '\n';
  //cout << "Config  : " << stringConfigFileName << "!" << '\n';
  if (isAcousticUsed) {
    implAcoustic.reset(new precice::SolverInterface(stringAccessorNameAcoustic,
                                                    stringConfigFileName,
                                                    *solverProcessIndex, *solverProcessSize));
  }
  if (isFluidUsed) {
    implFluid.reset(new precice::SolverInterface(stringAccessorNameFluid,
                                                 stringConfigFileName,
                                                 *solverProcessIndex, *solverProcessSize));
  }
  PRECICE_CHECK(implAcoustic != nullptr || implFluid != nullptr, "Either the Fluid interface or the Acoustic"
                                                                 " interface or both need to be used");
}

void precice_fastest_initialize_(
    double *   timestepLengthLimit,
    const int *useFluid)
{
  precice::impl::checkCorrectUsage(*useFluid);

  if (*useFluid == 0) {
    *timestepLengthLimit = implAcoustic->initialize();
  } else {
    *timestepLengthLimit = implFluid->initialize();
  }
}

void precice_fastest_initialize_data_(
    const int *useFluid)
{
  precice::impl::checkCorrectUsage(*useFluid);

  if (*useFluid == 0) {
    implAcoustic->initializeData();
  } else {
    implFluid->initializeData();
  }
}

void precice_fastest_advance_(
    double *   timestepLengthLimit,
    const int *useFluid)
{
  precice::impl::checkCorrectUsage(*useFluid);

  if (*useFluid == 0) {
    *timestepLengthLimit = implAcoustic->advance(*timestepLengthLimit);
  } else {
    *timestepLengthLimit = implFluid->advance(*timestepLengthLimit);
  }
}

void precice_fastest_finalize_(
    const int *useFluid)
{
  precice::impl::checkCorrectUsage(*useFluid);

  if (*useFluid == 0) {
    implAcoustic->finalize();
    implAcoustic.reset();
  } else {
    implFluid->finalize();
    implFluid.reset();
  }
}

void precice_fastest_action_required_(
    const char *action,
    int *       isRequired,
    const int * useFluid,
    int         lengthAction)
{
  precice::impl::checkCorrectUsage(*useFluid);

  int    strippedLength = precice::impl::strippedLength(action, lengthAction);
  string stringAction(action, strippedLength);
  if (*useFluid == 0) {
    if (implAcoustic->isActionRequired(stringAction)) {
      *isRequired = 1;
    } else {
      *isRequired = 0;
    }
  } else {
    if (implFluid->isActionRequired(stringAction)) {
      *isRequired = 1;
    } else {
      *isRequired = 0;
    }
  }
}

void precice_fastest_mark_action_fulfilled_(
    const char *action,
    const int * useFluid,
    int         lengthAction)
{
  precice::impl::checkCorrectUsage(*useFluid);

  int    strippedLength = precice::impl::strippedLength(action, lengthAction);
  string stringAction(action, strippedLength);
  if (*useFluid == 0) {
    implAcoustic->markActionFulfilled(stringAction);
  } else {
    implFluid->markActionFulfilled(stringAction);
  }
}

void precice_fastest_get_mesh_id_(
    const char *meshName,
    int *       meshID,
    const int * useFluid,
    int         lengthMeshName)
{
  precice::impl::checkCorrectUsage(*useFluid);

  int    strippedLength = precice::impl::strippedLength(meshName, lengthMeshName);
  string stringMeshName(meshName, strippedLength);
  if (*useFluid == 0) {
    *meshID = implAcoustic->getMeshID(stringMeshName);
  } else {
    *meshID = implFluid->getMeshID(stringMeshName);
  }
}

void precice_fastest_get_data_id_(
    const char *dataName,
    const int * meshID,
    int *       dataID,
    const int * useFluid,
    int         lengthDataName)
{
  precice::impl::checkCorrectUsage(*useFluid);

  int    strippedLength = precice::impl::strippedLength(dataName, lengthDataName);
  string stringDataName(dataName, strippedLength);
  if (*useFluid == 0) {
    *dataID = implAcoustic->getDataID(stringDataName, *meshID);
  } else {
    *dataID = implFluid->getDataID(stringDataName, *meshID);
  }
}

void precice_fastest_set_vertex_(
    const int *   meshID,
    const double *position,
    int *         vertexID,
    const int *   useFluid)
{
  precice::impl::checkCorrectUsage(*useFluid);

  if (*useFluid == 0) {
    *vertexID = implAcoustic->setMeshVertex(*meshID, position);
  } else {
    *vertexID = implFluid->setMeshVertex(*meshID, position);
  }
}

void precice_fastest_set_vertices_(
    const int *meshID,
    const int *size,
    double *   positions,
    int *      positionIDs,
    const int *useFluid)
{
  precice::impl::checkCorrectUsage(*useFluid);

  if (*useFluid == 0) {
    implAcoustic->setMeshVertices(*meshID, *size, positions, positionIDs);
  } else {
    implFluid->setMeshVertices(*meshID, *size, positions, positionIDs);
  }
}

void precice_fastest_write_bvdata_(
    const int *dataID,
    const int *size,
    int *      valueIndices,
    double *   values,
    const int *useFluid)
{
  precice::impl::checkCorrectUsage(*useFluid);

  if (*useFluid == 0) {
    implAcoustic->writeBlockVectorData(*dataID, *size, valueIndices, values);
  } else {
    implFluid->writeBlockVectorData(*dataID, *size, valueIndices, values);
  }
}

void precice_fastest_write_vdata_(
    const int *   dataID,
    const int *   valueIndex,
    const double *dataValue,
    const int *   useFluid)
{
  precice::impl::checkCorrectUsage(*useFluid);

  if (*useFluid == 0) {
    implAcoustic->writeVectorData(*dataID, *valueIndex, dataValue);
  } else {
    implFluid->writeVectorData(*dataID, *valueIndex, dataValue);
  }
}

void precice_fastest_write_bsdata_(
    const int *dataID,
    const int *size,
    int *      valueIndices,
    double *   values,
    const int *useFluid)
{
  precice::impl::checkCorrectUsage(*useFluid);

  if (*useFluid == 0) {
    implAcoustic->writeBlockScalarData(*dataID, *size, valueIndices, values);
  } else {
    implFluid->writeBlockScalarData(*dataID, *size, valueIndices, values);
  }
}

void precice_fastest_write_sdata_(
    const int *   dataID,
    const int *   valueIndex,
    const double *dataValue,
    const int *   useFluid)
{
  precice::impl::checkCorrectUsage(*useFluid);

  if (*useFluid == 0) {
    implAcoustic->writeScalarData(*dataID, *valueIndex, *dataValue);
  } else {
    implFluid->writeScalarData(*dataID, *valueIndex, *dataValue);
  }
}

void precice_fastest_read_bvdata_(
    const int *dataID,
    const int *size,
    int *      valueIndices,
    double *   values,
    const int *useFluid)
{
  precice::impl::checkCorrectUsage(*useFluid);

  if (*useFluid == 0) {
    implAcoustic->readBlockVectorData(*dataID, *size, valueIndices, values);
  } else {
    implFluid->readBlockVectorData(*dataID, *size, valueIndices, values);
  }
}

void precice_fastest_read_vdata_(
    const int *dataID,
    const int *valueIndex,
    double *   dataValue,
    const int *useFluid)
{
  precice::impl::checkCorrectUsage(*useFluid);

  if (*useFluid == 0) {
    implAcoustic->readVectorData(*dataID, *valueIndex, dataValue);
  } else {
    implFluid->readVectorData(*dataID, *valueIndex, dataValue);
  }
}

void precice_fastest_read_bsdata_(
    const int *dataID,
    const int *size,
    int *      valueIndices,
    double *   values,
    const int *useFluid)
{
  precice::impl::checkCorrectUsage(*useFluid);

  if (*useFluid == 0) {
    implAcoustic->readBlockScalarData(*dataID, *size, valueIndices, values);
  } else {
    implFluid->readBlockScalarData(*dataID, *size, valueIndices, values);
  }
}

void precice_fastest_read_sdata_(
    const int *dataID,
    const int *valueIndex,
    double *   dataValue,
    const int *useFluid)
{
  precice::impl::checkCorrectUsage(*useFluid);

  if (*useFluid == 0) {
    implAcoustic->readScalarData(*dataID, *valueIndex, *dataValue);
  } else {
    implFluid->readScalarData(*dataID, *valueIndex, *dataValue);
  }
}

void precice::impl::checkCorrectUsage(int useFluid)
{
  PRECICE_CHECK(useFluid == 0 || useFluid == 1, "useFluid needs to be either 0 or 1.");

  if (useFluid == 1) {
    PRECICE_CHECK(implFluid != nullptr, "The fluid interface has not been created properly. Be sure to call "
                                        "\"precicef_create\" with \"useFluid=1\" before any other call to preCICE if \"isFluid=1\".");
  } else if (useFluid == 0) {
    PRECICE_CHECK(implAcoustic != nullptr, "The acoustic interface has not been created properly. Be sure to call "
                                           "\"precicef_create\" with \"useAcoustic=1\" before any other call to preCICE if \"isFluid=0\".");
  }
}

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif
#ifdef __clang__
#pragma clang diagnostic pop
#endif
