#include <optional>
#include <string_view>

#include "cplscheme/Constants.hpp"
#include "precice/Participant.hpp"
#include "precice/impl/ParticipantImpl.hpp"
#include "precice/impl/versions.hpp"

namespace precice {

namespace {
std::string_view toSV(precice::string_view sv)
{
  // The given string_view may contain null chars.
  // We trim to the first null char here.
  std::string_view s{sv.data(), sv.size()};
  auto             trim_pos = s.find('\0');
  if (trim_pos != s.npos) {
    s.remove_suffix(s.size() - trim_pos);
  }
  return s;
}
} // namespace

/**
 * @brief Creates a Participant for the solver with the given name.
 *
 * @param participantName Name of the participant using the interface.
 * @param configurationFileName Name of the XML configuration file.
 * @param solverProcessIndex If the solver code runs with several processes, each process using preCICE has to specify its index, which has to start from 0 and end with solverProcessSize - 1.
 * @param solverProcessSize The number of solver processes using preCICE.
 *
 * \since 3.0.0
 * \note Renamed from SolverInterface in version 3.0.0
 */
Participant::Participant(
    ::precice::string_view participantName,
    ::precice::string_view configurationFileName,
    int                    solverProcessIndex,
    int                    solverProcessSize)
    : _impl(new impl::ParticipantImpl(toSV(participantName), toSV(configurationFileName), solverProcessIndex, solverProcessSize, std::nullopt))
{
}

/**
 * @brief Creates a Participant for the solver with the given name and a custom MPI communicator.
 *
 * @param participantName Name of the participant using the interface.
 * @param configurationFileName Name of the XML configuration file.
 * @param solverProcessIndex If the solver code runs with several processes, each process using preCICE has to specify its index, which has to start from 0 and end with solverProcessSize - 1.
 * @param solverProcessSize The number of solver processes using preCICE.
 * @param communicator A pointer to a custom MPI communicator.
 *
 * \since 1.6.0
 * \changed{3.0.0} Renamed from SolverInterface to Participant.
 */
Participant::Participant(
    ::precice::string_view participantName,
    ::precice::string_view configurationFileName,
    int                    solverProcessIndex,
    int                    solverProcessSize,
    void                  *communicator)
    : _impl(new impl::ParticipantImpl(toSV(participantName), toSV(configurationFileName), solverProcessIndex, solverProcessSize, {communicator}))
{
}

Participant::~Participant() = default;

/**
 * @brief Initializes preCICE and coupling data.
 *
 * \since 3.0.0
 * \note In older versions, this also called initializeData(). That function was removed in 3.0.0.
 */
void Participant::initialize()
{
  _impl->initialize();
}

/**
 * @brief Advances preCICE after the solver has computed one time step.
 *
 * @param computedTimeStepSize Size of time step used by the solver.
 *
 * \since 3.0.0
 */
void Participant::advance(
    double computedTimeStepSize)
{
  _impl->advance(computedTimeStepSize);
}

/**
 * @brief Finalizes preCICE.
 *
 * \since 3.0.0
 */
void Participant::finalize()
{
  return _impl->finalize();
}

/**
 * @brief Returns the spatial dimensionality of the given mesh.
 *
 * @param meshName the name of the mesh
 * @returns the dimensions of the given mesh
 *
 * \since 3.0.0
 * \note Replaces getDimensions() from version 2.x
 */
int Participant::getMeshDimensions(::precice::string_view meshName) const
{
  return _impl->getMeshDimensions(toSV(meshName));
}

/**
 * @brief Returns the spatial dimensionality of the given data on the given mesh.
 *
 * @param meshName the name of the mesh
 * @param dataName the name of the data
 * @returns the dimensions of the given data
 *
 * \since 3.0.0
 * \note Replaces getDimensions() from version 2.x
 */
int Participant::getDataDimensions(::precice::string_view meshName, ::precice::string_view dataName) const
{
  return _impl->getDataDimensions(toSV(meshName), toSV(dataName));
}

/**
 * @brief Checks if the coupled simulation is still ongoing.
 *
 * @returns whether the coupling is ongoing.
 *
 * \since 3.0.0
 */
bool Participant::isCouplingOngoing() const
{
  return _impl->isCouplingOngoing();
}

/**
 * @brief Checks if the current time window has been completed.
 *
 * @returns whether the current time window is complete.
 *
 * \since 2.0.0
 * \note Renamed from isTimestepComplete() in version 2.0.0.
 */
bool Participant::isTimeWindowComplete() const
{
  return _impl->isTimeWindowComplete();
}

/**
 * @brief Returns the maximum allowed time step size.
 *
 * @returns the maximum time step size.
 *
 * \since 3.0.0
 * \changed{3.0.0} Now returns 0.0 after the final advance(dt).
 */
double Participant::getMaxTimeStepSize() const
{
  return _impl->getMaxTimeStepSize();
}

/**
 * @brief Checks if initial data must be written before calling initialize().
 *
 * @returns whether initial data is required.
 *
 * \since 3.0.0
 * \note Replaces the old isActionRequired(actionWriteInitialData()) pattern.
 */
bool Participant::requiresInitialData()
{
  return _impl->requiresInitialData();
}

/**
 * @brief Checks if the participant must read a checkpoint.
 *
 * @returns whether reading a checkpoint is required.
 *
 * \since 3.0.0
 * \note Replaces the old isActionRequired(actionReadIterationCheckpoint()) pattern.
 */
bool Participant::requiresReadingCheckpoint()
{
  return _impl->requiresReadingCheckpoint();
}

/**
 * @brief Checks if the participant must write a checkpoint.
 *
 * @returns whether writing a checkpoint is required.
 *
 * \since 3.0.0
 * \note Replaces the old isActionRequired(actionWriteIterationCheckpoint()) pattern.
 */
bool Participant::requiresWritingCheckpoint()
{
  return _impl->requiresWritingCheckpoint();
}

/**
 * @brief Checks if the given mesh requires connectivity.
 *
 * @param meshName the name of the mesh
 * @returns whether mesh connectivity is required.
 *
 * \since 3.0.0
 * \note Renamed from isMeshConnectivityRequired() in version 3.0.0.
 */
bool Participant::requiresMeshConnectivityFor(::precice::string_view meshName) const
{
  return _impl->requiresMeshConnectivityFor(toSV(meshName));
}

/**
 * @brief Resets the mesh with the given name.
 *
 * @param meshName the name of the mesh to reset.
 *
 * \since 3.2.0
 * \note Part of experimental remeshing support added in 3.2.0.
 */
void Participant::resetMesh(::precice::string_view meshName)
{
  return _impl->resetMesh(toSV(meshName));
}

/**
 * @brief Checks if the given data on the given mesh requires gradient data.
 *
 * @param meshName the name of the mesh
 * @param dataName the name of the data
 * @returns whether gradient data is required.
 *
 * \since 3.0.0
 * \note Renamed from isGradientDataRequired() in version 3.0.0.
 */
bool Participant::requiresGradientDataFor(::precice::string_view meshName,
                                          ::precice::string_view dataName) const
{
  return _impl->requiresGradientDataFor(toSV(meshName), toSV(dataName));
}

/**
 * @brief Creates a mesh vertex and returns its ID.
 *
 * @param meshName the name of the mesh to add the vertex to.
 * @param coordinates a span of coordinates of the vertex.
 * @returns the ID of the created vertex.
 *
 * \since 3.0.0
 * \changed{3.0.0} Now uses precice::span instead of raw pointers.
 */
VertexID Participant::setMeshVertex(
    ::precice::string_view        meshName,
    ::precice::span<const double> coordinates)
{
  return _impl->setMeshVertex(toSV(meshName), coordinates);
}

/**
 * @brief Returns the number of vertices of a mesh.
 *
 * @param meshName the name of the mesh.
 * @returns the number of vertices of the mesh.
 *
 * \since 3.0.0
 */
int Participant::getMeshVertexSize(
    ::precice::string_view meshName) const
{
  return _impl->getMeshVertexSize(toSV(meshName));
}

/**
 * @brief Creates multiple mesh vertices.
 *
 * @param meshName the name of the mesh to add the vertices to.
 * @param coordinates a span of coordinates of the vertices.
 * @param ids a span to write the IDs of the created vertices to.
 *
 * \since 3.0.0
 * \changed{3.0.0} Now uses precice::span instead of raw pointers.
 */
void Participant::setMeshVertices(
    ::precice::string_view        meshName,
    ::precice::span<const double> coordinates,
    ::precice::span<VertexID>     ids)
{
  _impl->setMeshVertices(toSV(meshName), coordinates, ids);
}

/**
 * @brief Sets a mesh edge between two vertices.
 *
 * @param meshName the name of the mesh.
 * @param first the ID of the first vertex.
 * @param second the ID of the second vertex.
 *
 * \since 3.0.0
 * \changed{3.0.0} No longer returns an edge ID.
 */
void Participant::setMeshEdge(
    ::precice::string_view meshName,
    VertexID               first,
    VertexID               second)
{
  _impl->setMeshEdge(toSV(meshName), first, second);
}

/**
 * @brief Sets multiple mesh edges at once.
 *
 * @param meshName the name of the mesh.
 * @param ids a span of vertex IDs in pairs (first, second).
 *
 * \since 3.0.0
 */
void Participant::setMeshEdges(
    ::precice::string_view          meshName,
    ::precice::span<const VertexID> ids)
{
  _impl->setMeshEdges(toSV(meshName), ids);
}

/**
 * @brief Sets a triangle defined by three vertex IDs.
 *
 * @param meshName the name of the mesh.
 * @param first the ID of the first vertex.
 * @param second the ID of the second vertex.
 * @param third the ID of the third vertex.
 *
 * \since 3.0.0
 * \changed{3.0.0} Now accepts vertex IDs instead of edge IDs.
 */
void Participant::setMeshTriangle(
    ::precice::string_view meshName,
    VertexID               first,
    VertexID               second,
    VertexID               third)
{
  _impl->setMeshTriangle(toSV(meshName), first, second, third);
}

/**
 * @brief Sets multiple mesh triangles at once.
 *
 * @param meshName the name of the mesh.
 * @param ids a span of vertex IDs in triplets.
 *
 * \since 3.0.0
 */
void Participant::setMeshTriangles(
    ::precice::string_view          meshName,
    ::precice::span<const VertexID> ids)
{
  _impl->setMeshTriangles(toSV(meshName), ids);
}

/**
 * @brief Sets a quad defined by four vertex IDs.
 *
 * @param meshName the name of the mesh.
 * @param first the ID of the first vertex.
 * @param second the ID of the second vertex.
 * @param third the ID of the third vertex.
 * @param fourth the ID of the fourth vertex.
 *
 * \since 3.0.0
 * \changed{3.3.0} Added support for 2D meshes.
 */
void Participant::setMeshQuad(
    ::precice::string_view meshName,
    VertexID               first,
    VertexID               second,
    VertexID               third,
    VertexID               fourth)
{
  _impl->setMeshQuad(toSV(meshName), first, second, third, fourth);
}

/**
 * @brief Sets multiple mesh quads at once.
 *
 * @param meshName the name of the mesh.
 * @param ids a span of vertex IDs in groups of four.
 *
 * \since 3.0.0
 * \changed{3.3.0} Added support for 2D meshes.
 */
void Participant::setMeshQuads(
    ::precice::string_view          meshName,
    ::precice::span<const VertexID> ids)
{
  _impl->setMeshQuads(toSV(meshName), ids);
}

/**
 * @brief Sets a tetrahedron defined by four vertex IDs.
 *
 * @param meshName the name of the mesh.
 * @param first the ID of the first vertex.
 * @param second the ID of the second vertex.
 * @param third the ID of the third vertex.
 * @param fourth the ID of the fourth vertex.
 *
 * \since 3.0.0
 */
void Participant::setMeshTetrahedron(
    ::precice::string_view meshName,
    VertexID               first,
    VertexID               second,
    VertexID               third,
    VertexID               fourth)
{
  _impl->setMeshTetrahedron(toSV(meshName), first, second, third, fourth);
}

/**
 * @brief Sets multiple mesh tetrahedra at once.
 *
 * @param meshName the name of the mesh.
 * @param ids a span of vertex IDs in groups of four.
 *
 * \since 3.0.0
 */
void Participant::setMeshTetrahedra(
    ::precice::string_view          meshName,
    ::precice::span<const VertexID> ids)
{
  _impl->setMeshTetrahedra(toSV(meshName), ids);
}

/**
 * @brief Writes data to a mesh.
 *
 * @param meshName the name of the mesh.
 * @param dataName the name of the data.
 * @param ids a span of vertex IDs to write data to.
 * @param values a span of values to write.
 *
 * \since 3.0.0
 * \changed{3.0.0} Replaces writeBlockScalarData and writeBlockVectorData. Now uses precice::span.
 */
void Participant::writeData(
    ::precice::string_view          meshName,
    ::precice::string_view          dataName,
    ::precice::span<const VertexID> ids,
    ::precice::span<const double>   values)
{
  _impl->writeData(toSV(meshName), toSV(dataName), ids, values);
}

/**
 * @brief Reads data from a mesh.
 *
 * @param meshName the name of the mesh.
 * @param dataName the name of the data.
 * @param ids a span of vertex IDs to read data from.
 * @param relativeReadTime relative time to read data at. Use 0.0 for current time, getMaxTimeStepSize() for end of time window.
 * @param values a span to write the read values to.
 *
 * \since 3.0.0
 * \changed{3.0.0} relativeReadTime is now mandatory. Replaces readBlockScalarData and readBlockVectorData.
 */
void Participant::readData(
    ::precice::string_view          meshName,
    ::precice::string_view          dataName,
    ::precice::span<const VertexID> ids,
    double                          relativeReadTime,
    ::precice::span<double>         values) const
{
  _impl->readData(toSV(meshName), toSV(dataName), ids, relativeReadTime, values);
}

/**
 * @brief Maps data and reads it for the given coordinates.
 *
 * Experimental API for just-in-time (JIT) mapping. Allows reading data at arbitrary coordinates
 * without defining a mesh beforehand.
 *
 * @param meshName the name of the mesh.
 * @param dataName the name of the data.
 * @param coordinates a span of coordinates to read data at.
 * @param relativeReadTime relative time to read data at.
 * @param values a span to write the read values to.
 *
 * \since 3.2.0
 */
void Participant::mapAndReadData(
    ::precice::string_view        meshName,
    ::precice::string_view        dataName,
    ::precice::span<const double> coordinates,
    double                        relativeReadTime,
    ::precice::span<double>       values) const
{
  _impl->mapAndReadData(toSV(meshName), toSV(dataName), coordinates, relativeReadTime, values);
}

/**
 * @brief Writes data and maps it for the given coordinates.
 *
 * Experimental API for just-in-time (JIT) mapping. Allows writing data at arbitrary coordinates
 * without defining a mesh beforehand.
 *
 * @param meshName the name of the mesh.
 * @param dataName the name of the data.
 * @param coordinates a span of coordinates to write data at.
 * @param values a span of values to write.
 *
 * \since 3.2.0
 */
void Participant::writeAndMapData(
    ::precice::string_view        meshName,
    ::precice::string_view        dataName,
    ::precice::span<const double> coordinates,
    ::precice::span<const double> values)
{
  _impl->writeAndMapData(toSV(meshName), toSV(dataName), coordinates, values);
}

/**
 * @brief Sets a bounding box for the mesh access region.
 *
 * Required when using direct mesh access (api-access). Must be called before initialize().
 *
 * @param meshName the name of the mesh.
 * @param boundingBox a span defining the bounding box.
 *
 * \since 3.0.0
 * \changed{3.2.0} Now requires api-access="true" in the configuration. Can now be called once per mesh instead of once per participant.
 */
void Participant::setMeshAccessRegion(::precice::string_view        meshName,
                                      ::precice::span<const double> boundingBox) const
{
  _impl->setMeshAccessRegion(toSV(meshName), boundingBox);
}

/**
 * @brief Gets vertex IDs and coordinates of a mesh.
 *
 * Used together with setMeshAccessRegion() for direct mesh access.
 *
 * @param meshName the name of the mesh.
 * @param ids a span to write the vertex IDs to.
 * @param coordinates a span to write the vertex coordinates to.
 *
 * \since 3.0.0
 * \changed{3.0.0} Renamed from getMeshVerticesAndIDs().
 * \changed{3.2.0} Now requires api-access="true" in the configuration. Returns only vertices within the defined access region.
 */
void Participant::getMeshVertexIDsAndCoordinates(::precice::string_view    meshName,
                                                 ::precice::span<VertexID> ids,
                                                 ::precice::span<double>   coordinates) const
{
  _impl->getMeshVertexIDsAndCoordinates(toSV(meshName), ids, coordinates);
}

/**
 * @brief Writes gradient data to a mesh.
 *
 * @param meshName the name of the mesh.
 * @param dataName the name of the data.
 * @param ids a span of vertex IDs to write gradient data to.
 * @param gradients a span of gradient values to write.
 *
 * \since 3.0.0
 * \note Gradient data support was added as experimental in 2.4.0 and finalized in 3.0.0.
 */
void Participant::writeGradientData(
    ::precice::string_view          meshName,
    ::precice::string_view          dataName,
    ::precice::span<const VertexID> ids,
    ::precice::span<const double>   gradients)
{
  _impl->writeGradientData(toSV(meshName), toSV(dataName), ids, gradients);
}

/**
 * @brief Starts a user-defined profiling section.
 *
 * Allows measuring and displaying profiling of adapter code in the context of the whole participant or simulation.
 *
 * @param sectionName the name of the profiling section.
 *
 * \since 3.2.0
 */
void Participant::startProfilingSection(::precice::string_view sectionName)
{
  _impl->startProfilingSection(toSV(sectionName));
}

/**
 * @brief Stops the last started user-defined profiling section.
 *
 * \since 3.2.0
 */
void Participant::stopLastProfilingSection()
{
  _impl->stopLastProfilingSection();
}

} // namespace precice
