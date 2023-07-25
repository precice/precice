#pragma once

#include <memory>
#include <precice/Version.h>
#include <precice/export.h>
#include <precice/span.hpp>
#include "precice/types.hpp"

/**
 * forward declarations.
 */
namespace precice {
namespace impl {
class ParticipantImpl;
}
namespace testing {
struct WhiteboxAccessor;
}
} // namespace precice

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {

/// convenience type for compatibility
using string_view = ::precice::span<const char>;

/**
 * @brief Main Application Programming Interface of preCICE
 *
 * To adapt a solver to preCICE, follow the following main structure:
 *
 * -# Create an object of Participant with Participant()
 * -# Initialize preCICE with Participant::initialize()
 * -# Advance to the next (time)step with Participant::advance()
 * -# Finalize preCICE with Participant::finalize()
 *
 *  @note
 *  We use solver, simulation code, and participant as synonyms.
 *  The preferred name in the documentation is participant.
 */
class PRECICE_API Participant {
public:
  ///@name Construction and Configuration
  ///@{

  /**
   * @brief Constructs a Participant for the given participant
   *
   * @param[in] participantName Name of the participant using the interface. Has to
   *        match the name given for a participant in the xml configuration file.
   * @param[in] configurationFileName Name (with path) of the xml configuration file.
   * @param[in] solverProcessIndex If the solver code runs with several processes,
   *        each process using preCICE has to specify its index, which has to start
   *        from 0 and end with solverProcessSize - 1.
   * @param[in] solverProcessSize The number of solver processes using preCICE.
   */
  Participant(
      ::precice::string_view participantName,
      ::precice::string_view configurationFileName,
      int                    solverProcessIndex,
      int                    solverProcessSize);

  /**
   * @brief Constructs a Participant for the given participant and a custom MPI communicator.
   *
   * @param[in] participantName Name of the participant using the interface. Has to
   *        match the name given for a participant in the xml configuration file.
   * @param[in] configurationFileName Name (with path) of the xml configuration file.
   * @param[in] solverProcessIndex If the solver code runs with several processes,
   *        each process using preCICE has to specify its index, which has to start
   *        from 0 and end with solverProcessSize - 1.
   * @param[in] solverProcessSize The number of solver processes using preCICE.
   * @param[in] communicator A pointer to an MPI_Comm to use as communicator.
   */
  Participant(
      ::precice::string_view participantName,
      ::precice::string_view configurationFileName,
      int                    solverProcessIndex,
      int                    solverProcessSize,
      void *                 communicator);

  ~Participant();

  ///@}

  /// @name Steering Methods
  ///@{

  /**
   * @brief Fully initializes preCICE and coupling data.
   *
   * - Sets up a connection to the other participants of the coupled simulation.
   * - Creates all meshes, solver meshes need to be submitted before.
   * - Receives first coupling data. The starting values for coupling data are zero by default.
   * - Determines maximum allowed size of the first time step to be computed.
   *
   * @see getMaxTimeStepSize
   *
   * @pre initialize() has not yet been called.
   *
   * @post Parallel communication to the coupling partner(s) is setup.
   * @post Meshes are exchanged between coupling partners and the parallel partitions are created.
   * @post Initial coupling data was exchanged.
   */
  void initialize();

  /**
   * @brief Advances preCICE after the solver has computed one time step.
   *
   * - Sends and resets coupling data written by solver to coupling partners.
   * - Receives coupling data read by solver.
   * - Computes and applies data mappings.
   * - Computes acceleration of coupling data.
   * - Exchanges and computes information regarding the state of the coupled
   *   simulation.
   *
   * @param[in] computedTimeStepSize Size of time step used by the solver.
   *
   * @see getMaxTimeStepSize to get the maximum allowed value for \p computedTimeStepSize.
   *
   * @pre initialize() has been called successfully.
   * @pre The solver has computed one time step.
   * @pre The solver has written all coupling data.
   * @pre isCouplngOngoing() returns true.
   * @pre finalize() has not yet been called.
   *
   * @post Coupling data values specified in the configuration are exchanged.
   * @post Coupling scheme state (computed time, computed time steps, ...) is updated.
   * @post The coupling state is logged.
   * @post Configured data mapping schemes are applied.
   * @post [Second Participant] Configured acceleration schemes are applied.
   * @post Meshes with data are exported to files if configured.
   */
  void advance(double computedTimeStepSize);

  /**
   * @brief Finalizes preCICE.
   *
   * If initialize() has been called:
   *
   * - Synchronizes with remote participants
   * - handles final exports
   * - cleans up general state
   *
   * Always:
   *
   * - flushes and finalizes Events
   * - finalizes managed PETSc
   * - finalizes managed MPI
   *
   * @pre finalize() has not been called.
   *
   * @post Communication channels are closed.
   * @post Meshes and data are deallocated
   * @post Finalized managed PETSc
   * @post Finalized managed MPI
   *
   * @see isCouplingOngoing()
   */
  void finalize();

  ///@}

  ///@name Status Queries
  ///@{

  /**
   * @brief Returns the spatial dimensionality of the given mesh.
   *
   * @param[in] meshName the name of the associated mesh
   *
   * @returns the dimensions of the given mesh
   */
  int getMeshDimensions(::precice::string_view meshName) const;

  /**
   * @brief Returns the spatial dimensionality of the given data on the given mesh.
   *
   * Note that vectorial data dimensionality directly depends on the spacial dimensionality of the mesh.
   *
   * @param[in] meshName the name of the associated mesh
   * @param[in] dataName the name of the data to get the dimensions for
   *
   * @returns the dimensions of the given Data
   *
   * @see getMeshDimensions
   */
  int getDataDimensions(::precice::string_view meshName, ::precice::string_view dataName) const;

  /**
   * @brief Checks if the coupled simulation is still ongoing.
   *
   * @returns whether the coupling is ongoing.
   *
   * A coupling is ongoing as long as
   * - the maximum number of time windows has not been reached, and
   * - the final time has not been reached.
   *
   * @pre initialize() has been called successfully.
   *
   * @see advance()
   *
   * @note
   * The user should call finalize() after this function returns false.
   */
  bool isCouplingOngoing() const;

  /**
   * @brief Checks if the current coupling window is completed.
   *
   * @returns whether the current time window is complete.
   *
   * The following reasons require several solver time steps per time window:
   * - A solver chooses to perform subcycling, i.e. using a smaller timestep
   *   than the time window.
   * - An implicit coupling iteration is not yet converged.
   *
   * Hence, a time window is complete if we reach the end of the time window
   * and the implicit coupling has converged.
   *
   * For implicit coupling this condition is equivalent with the requirement to
   * write an iteration checkpoint. This is, however, not the case for explicit
   * coupling.
   *
   * @pre initialize() has been called successfully.
   */
  bool isTimeWindowComplete() const;

  /**
   * @brief Get the maximum allowed time step size of the current window.
   *
   * @returns Maximum size of time step to be computed by solver.
   *
   * Allows the user to query the maximum allowed time step size in the current window.
   * This should be used to compute the actual time step that the solver uses.
   *
   * @pre initialize() has been called successfully.
   */
  double getMaxTimeStepSize() const;

  ///@}

  ///@name Requirements
  ///@{

  /** Checks if the participant is required to provide initial data.
   *
   * If true, then the participant needs to write initial data to defined vertices
   * prior to calling initialize().
   *
   * @pre initialize() has not yet been called
   */
  bool requiresInitialData();

  /** Checks if the participant is required to write an iteration checkpoint.
   *
   * If true, the participant is required to write an iteration checkpoint before
   * calling advance().
   *
   * preCICE refuses to proceed if writing a checkpoint is required,
   * but this method isn't called prior to advance().
   *
   * @pre initialize() has been called
   *
   * @see requiresReadingCheckpoint()
   */
  bool requiresWritingCheckpoint();

  /** Checks if the participant is required to read an iteration checkpoint.
   *
   * If true, the participant is required to read an iteration checkpoint before
   * calling advance().
   *
   * preCICE refuses to proceed if reading a checkpoint is required,
   * but this method isn't called prior to advance().
   *
   * @note This function returns false before the first call to advance().
   *
   * @pre initialize() has been called
   *
   * @see requiresWritingCheckpoint()
   */
  bool requiresReadingCheckpoint();

  ///@}

  /** @name Mesh Access
   * @anchor precice-mesh-access
   *
   * Connectivity is optional.
   * Use requiresMeshConnectivityFor() to check if the current participant can make use of the connectivity.
   *
   * Only set the mesh connectivity that you require for your use-case and use the face/cell elements that your solver provides.
   * preCICE removes all connectivity duplicates in initialize().
   *
   * We recommend you to do the following depending on your case:
   *
   * - **2D surface coupling:** Use setMeshEdge() and setMeshEdges() to specify the coupling interface.
   * - **2D volume coupling:** Use setMeshTriangle() and setMeshTriangles() to specify the coupling area.
   * - **3D surface coupling:** Use setMeshTriangle() and setMeshTriangles() to specify the coupling interface.
   * - **3D volume coupling:** Use setMeshTetrahedron() and setMeshTetrahedra() to specify the coupling volume.
   *
   * As an alternative to triangles, preCICE supports **planar** quads using setMeshQuad() and setMeshQuads().
   * These quads will be triangulated by preCICE, hence specifying triangles is generally preferred.
   * Before using quads, we recommended to check if your solver provides a way to traverse triangulated faces.
   *
   *@{
   */

  /*
   * @brief Resets mesh with given ID.
   *
   * @experimental
   *
   * Has to be called, every time the positions for data to be mapped
   * changes. Only has an effect, if the mapping used is non-stationary and
   * non-incremental.
   */
  //  void resetMesh ( ::precice::string_view meshName );

  /**
   * @brief Checks if the mesh with given name is used by a solver.
   *
   * @param[in] meshName the name of the mesh
   * @returns whether the mesh is used.
   */
  bool hasMesh(::precice::string_view meshName) const;

  /**
   * @brief Checks if the given mesh requires connectivity.
   *
   * preCICE may require connectivity information from the solver and
   * ignores any API calls regarding connectivity if it is not required.
   * Use this function to conditionally generate this connectivity.
   *
   * @param[in] meshName the name of the mesh
   * @returns whether connectivity is required
   */
  bool requiresMeshConnectivityFor(::precice::string_view meshName) const;

  /**
   * @brief Creates a mesh vertex
   *
   * @param[in] meshName the name of the mesh to add the vertex to.
   * @param[in] position the coordinates of the vertex.
   * @returns the id of the created vertex
   *
   * @pre initialize() has not yet been called
   * @pre position.size() == getMeshDimensions(meshName)
   *
   * @see getMeshDimensions()
   */
  int setMeshVertex(
      ::precice::string_view        meshName,
      ::precice::span<const double> position);

  /**
   * @brief Returns the number of vertices of a mesh.
   *
   * @param[in] meshName the name of the mesh
   * @returns the amount of the vertices of the mesh
   *
   * @pre This function can be called on received meshes as well as provided
   * meshes. However, you need to call this function after @p initialize(),
   * if the \p meshName corresponds to a received mesh, since the relevant mesh data
   * is exchanged during the @p initialize() call.
   */
  int getMeshVertexSize(::precice::string_view meshName) const;

  /**
   * @brief Creates multiple mesh vertices
   *
   * @param[in] meshName the name of the mesh to add the vertices to.
   * @param[in] positions a span to the coordinates of the vertices
   *            The 2D-format is (d0x, d0y, d1x, d1y, ..., dnx, dny)
   *            The 3D-format is (d0x, d0y, d0z, d1x, d1y, d1z, ..., dnx, dny, dnz)
   *
   * @param[out] ids The ids of the created vertices
   *
   * @pre initialize() has not yet been called
   * @pre position.size() == getMeshDimensions(meshName) * ids.size()
   *
   * @see getDimensions()
   */
  void setMeshVertices(
      ::precice::string_view        meshName,
      ::precice::span<const double> positions,
      ::precice::span<VertexID>     ids);

  /**
   * @brief Sets a mesh edge from vertex IDs
   *
   * Ignored if preCICE doesn't require connectivity for the mesh.
   *
   * @note The order of vertices does not matter.
   *
   * @param[in] meshName name of the mesh to add the edge to
   * @param[in] firstVertexID ID of the first vertex of the edge
   * @param[in] secondVertexID ID of the second vertex of the edge
   *
   * @pre vertices with firstVertexID and secondVertexID were added to the mesh with the name meshName
   */
  void setMeshEdge(
      ::precice::string_view meshName,
      int                    firstVertexID,
      int                    secondVertexID);

  /**
   * @brief Sets multiple mesh edge from vertex IDs
   *
   * vertices contain pairs of vertex indices for each edge to define.
   * The format follows: e1a, e1b, e2a, e2b, ...
   * Ignored if preCICE doesn't require connectivity for the mesh.
   *
   * @note The order of vertices per edge does not matter.
   *
   * @param[in] meshName the name of the mesh to add the edges to
   * @param[in] vertices an array containing 2*size vertex IDs
   *
   * @pre vertices were added to the mesh with the name meshName
   * @pre vertices.size() is multiple of 2
   *
   * @see requiresMeshConnectivityFor()
   */
  void setMeshEdges(
      ::precice::string_view          meshName,
      ::precice::span<const VertexID> vertices);

  /**
   * @brief Sets mesh triangle from vertex IDs.
   *
   *
   *
   * @note The order of vertices does not matter.
   *
   * @param[in] meshName name of the mesh to add the triangle to
   * @param[in] firstVertexID ID of the first vertex of the triangle
   * @param[in] secondVertexID ID of the second vertex of the triangle
   * @param[in] thirdVertexID ID of the third vertex of the triangle
   *
   * @pre edges with firstVertexID, secondVertexID, and thirdVertexID were added to the mesh with the name meshName
   *
   * @see requiresMeshConnectivityFor()
   */
  void setMeshTriangle(
      ::precice::string_view meshName,
      int                    firstVertexID,
      int                    secondVertexID,
      int                    thirdVertexID);

  /**
   * @brief Sets multiple mesh triangles from vertex IDs
   *
   * vertices contain triples of vertex indices for each triangle to define.
   * The format follows: t1a, t1b, t1c, t2a, t2b, t2c, ...
   * Ignored if preCICE doesn't require connectivity for the mesh.
   *
   * @note The order of vertices per triangle does not matter.
   *
   * @param[in] meshName name of the mesh to add the triangles to
   * @param[in] vertices an array containing 3*size vertex IDs
   *
   * @pre vertices were added to the mesh with the name meshName
   * @pre vertices.size() is multiple of 3
   *
   * @see requiresMeshConnectivityFor()
   */
  void setMeshTriangles(
      ::precice::string_view          meshName,
      ::precice::span<const VertexID> vertices);

  /**
   * @brief Sets a planar surface mesh quadrangle from vertex IDs.
   *
   * The planar quad will be triangulated, maximizing area-to-circumference.
   * Ignored if preCICE doesn't require connectivity for the mesh.
   *
   * @warning The order of vertices does not matter, however, only planar quads are allowed.
   *
   * @param[in] meshName name of the mesh to add the Quad to
   * @param[in] firstVertexID ID of the first vertex of the Quad
   * @param[in] secondVertexID ID of the second vertex of the Quad
   * @param[in] thirdVertexID ID of the third vertex of the Quad
   * @param[in] fourthVertexID ID of the fourth vertex of the Quad
   *
   * @pre vertices with firstVertexID, secondVertexID, thirdVertexID, and fourthVertexID were added to the mesh with the name meshName
   *
   * @see requiresMeshConnectivityFor()
   */
  void setMeshQuad(
      ::precice::string_view meshName,
      int                    firstVertexID,
      int                    secondVertexID,
      int                    thirdVertexID,
      int                    fourthVertexID);

  /**
   * @brief Sets multiple mesh quads from vertex IDs
   *
   * vertices contain quadruples of vertex indices for each quad to define.
   * The format follows: q1a, q1b, q1c, q1d, q2a, q2b, q2c, q2d, ...
   *
   * Each planar quad will be triangulated, maximizing area-to-circumference.
   * Ignored if preCICE doesn't require connectivity for the mesh.
   *
   * @warning The order of vertices per quad does not matter, however, only planar quads are allowed.
   *
   * @param[in] meshName name of the mesh to add the quad to
   * @param[in] vertices an array containing 4*size vertex IDs
   *
   * @pre vertices were added to the mesh with the name meshName
   * @pre vertices.size() is multiple of 4
   *
   * @see requiresMeshConnectivityFor()
   */
  void setMeshQuads(
      ::precice::string_view          meshName,
      ::precice::span<const VertexID> vertices);

  /**
   * @brief Set tetrahedron in 3D mesh from vertex ID
   *
   * @note The order of vertices does not matter.
   *
   * @param[in] meshName name of the mesh to add the Tetrahedron to
   * @param[in] firstVertexID ID of the first vertex of the Tetrahedron
   * @param[in] secondVertexID ID of the second vertex of the Tetrahedron
   * @param[in] thirdVertexID ID of the third vertex of the Tetrahedron
   * @param[in] fourthVertexID ID of the fourth vertex of the Tetrahedron
   *
   * @pre vertices with firstVertexID, secondVertexID, thirdVertexID, and fourthVertexID were added to the mesh with the name meshName
   *
   * @see requiresMeshConnectivityFor()
   */
  void setMeshTetrahedron(
      ::precice::string_view meshName,
      int                    firstVertexID,
      int                    secondVertexID,
      int                    thirdVertexID,
      int                    fourthVertexID);

  /**
   * @brief Sets multiple mesh tetrahedra from vertex IDs
   *
   * vertices contain quadruples of vertex indices for each tetrahedron to define.
   * The format follows: t1a, t1b, t1c, t1d, t2a, t2b, t2c, t2d, ...
   * Ignored if preCICE doesn't require connectivity for the mesh.
   *
   * @note The order of vertices per tetrahedron does not matter.
   *
   * @param[in] meshName name of the mesh to add the tetrahedra to
   * @param[in] vertices an array containing 4*size vertex IDs
   *
   * @pre vertices were added to the mesh with the name meshName
   * @pre vertices.size() is multiple of 4
   *
   * @see requiresMeshConnectivityFor()
   */
  void setMeshTetrahedra(
      ::precice::string_view          meshName,
      ::precice::span<const VertexID> vertices);

  ///@}

  ///@name Data Access
  ///@{

  /**
   * @brief Checks if the data with given name is used by a solver and mesh.
   *
   * @param[in] meshName the name of the associated mesh
   * @param[in] dataName the name of the data to check
   * @returns whether the mesh contains the data.
   */
  bool hasData(
      ::precice::string_view meshName,
      ::precice::string_view dataName) const;

  /**
   * @brief Writes data to a mesh.
   *
   * This function writes values of specified vertices to data of a mesh.
   * Values are provided as a block of continuous memory defined by values.
   * The order of the provided data follows the order specified by vertices.
   *
   * The 1D/Scalar-format of values is (d0, d1, ..., dn)
   * The 2D-format of values is (d0x, d0y, d1x, d1y, ..., dnx, dny)
   * The 3D-format of values is (d0x, d0y, d0z, d1x, d1y, d1z, ..., dnx, dny, dnz)
   *
   * @param[in] meshName the name of mesh that hold the data.
   * @param[in] dataName the name of the data to write to.
   * @param[in] vertices the vertex ids of the vertices to write data to.
   * @param[in] values the values to write to preCICE.
   *
   * @pre every VertexID in vertices is a return value of setMeshVertex or setMeshVertices
   * @pre values.size() == getDataDimensions(meshName, dataName) * vertices.size()
   *
   * @see Participant::setMeshVertex()
   * @see Participant::setMeshVertices()
   * @see Participant::getDataDimensions()
   */
  void writeData(
      ::precice::string_view          meshName,
      ::precice::string_view          dataName,
      ::precice::span<const VertexID> vertices,
      ::precice::span<const double>   values);

  /**
   * @brief Reads data values from a mesh. Values correspond to a given point in time relative to the beginning of the current timestep.
   *
   * This function reads values of specified vertices from data of a mesh.
   * Values are read into a block of continuous memory defined by values in the order specified by vertices.
   *
   * The 1D/Scalar-format of values is (d0, d1, ..., dn)
   * The 2D-format of values is (d0x, d0y, d1x, d1y, ..., dnx, dny)
   * The 3D-format of values is (d0x, d0y, d0z, d1x, d1y, d1z, ..., dnx, dny, dnz)
   *
   * The data is read at relativeReadTime, which indicates the point in time measured from the beginning of the current time step.
   * relativeReadTime = 0 corresponds to data at the beginning of the time step. Assuming that the user will call advance(dt) at the
   * end of the time step, dt indicates the size of the current time step. Then relativeReadTime = dt corresponds to the data at
   * the end of the time step.
   *
   * @param[in] meshName the name of mesh that hold the data.
   * @param[in] dataName the name of the data to read from.
   * @param[in] vertices the vertex ids of the vertices to read data from.
   * @param[in] relativeReadTime Point in time where data is read relative to the beginning of the current time step.
   * @param[out] values the destination memory to read the data from.
   *
   * @pre every VertexID in vertices is a return value of setMeshVertex or setMeshVertices
   * @pre values.size() == getDataDimensions(meshName, dataName) * vertices.size()
   *
   * @post values contain the read data as specified in the above format.
   *
   * @see Participant::setMeshVertex()
   * @see Participant::setMeshVertices()
   * @see Participant::getDataDimensions()
   */
  void readData(
      ::precice::string_view          meshName,
      ::precice::string_view          dataName,
      ::precice::span<const VertexID> vertices,
      double                          relativeReadTime,
      ::precice::span<double>         values) const;

  ///@}

  /** @name Experimental: Direct Access
   * These API functions are \b experimental and may change in future versions.
   */
  ///@{

  /**
   * @brief setMeshAccessRegion Define a region of interest on a received mesh
   *        (<receive-mesh ... from="otherParticipant />") in order to receive
   *        only a certain mesh region. Have a look at the website under
   *        https://precice.org/couple-your-code-direct-access.html or
   *        navigate manually to the page  Docs->Couple your code
   *        -> Advanced topics -> Accessing received meshes directly for
   *        a comprehensive documentation
   *
   * @experimental
   *
   * This function is required if you don't want to use the mapping
   * schemes in preCICE, but rather want to use your own solver for
   * data mapping. As opposed to the usual preCICE mapping, only a
   * single mesh (from the other participant) is now involved in this
   * situation since an 'own' mesh defined by the participant itself
   * is not required any more. In order to re-partition the received
   * mesh, the participant needs to define the mesh region it wants
   * read data from and write data to.
   * The mesh region is specified through an axis-aligned bounding
   * box given by the lower and upper [min and max] bounding-box
   * limits in each space dimension [x, y, z].
   *
   * @note Defining a bounding box for serial runs of the solver (not
   * to be confused with serial coupling mode) is valid. However, a
   * warning is raised in case vertices are filtered out completely
   * on the receiving side, since the associated data values of the
   * filtered vertices are filled with zero data.
   *
   * @note This function can only be called once per participant and
   * rank and trying to call it more than once results in an error.
   *
   * @note If you combine the direct access with a mpping (say you want
   * to read data from a defined mesh, as usual, but you want to directly
   * access and write data on a received mesh without a mapping) you may
   * not need this function at all since the region of interest is already
   * defined through the defined mesh used for data reading. This is the
   * case if you define any mapping involving the directly accessed mesh
   * on the receiving participant. (In parallel, only the cases
   * read-consistent and write-conservative are relevant, as usual).
   *
   * @note The safety factor scaling (see safety-factor in the configuration
   * file) is not applied to the defined access region and a specified safety
   * will be ignored in case there is no additional mapping involved. However,
   * in case a mapping is in addition to the direct access involved, you will
   * receive (and gain access to) vertices inside the defined access region
   * plus vertices inside the safety factor region resulting from the mapping.
   * The default value of the safety factor is 0.5,i.e., the defined access
   * region as computed through the involved provided mesh is by 50% enlarged.
   *
   * @param[in] meshName name of the mesh you want to access through the bounding box
   * @param[in] boundingBox Axis aligned bounding boxes which has in 3D the format
   *            [x_min, x_max, y_min, y_max, z_min, z_max]
   *
   * @pre @p initialize() has not yet been called.
   * @pre boundingBox.size() == 2 * getMeshDimensions(meshName)
   */
  void setMeshAccessRegion(
      ::precice::string_view        meshName,
      ::precice::span<const double> boundingBox) const;

  /**
   * @brief getMeshVertexIDsAndCoordinates Iterates over the region of
   *        interest defined by bounding boxes and reads the corresponding
   *        coordinates omitting the mapping.
   *
   * @experimental
   *
   * @param[in]  meshName corresponding mesh name
   * @param[out] ids ids corresponding to the coordinates
   * @param[out] coordinates the coordinates associated to the \p ids and
   *             corresponding data values
   *
   * @pre ids.size() == getMeshVertexSize(meshName)
   * @pre coordinates.size() == getMeshVertexSize(meshName) * getMeshDimensions(meshName)
   *
   * @pre This function can be called on received meshes as well as provided
   * meshes. However, you need to call this function after @p initialize(),
   * if the \p meshName corresponds to a received mesh, since the relevant mesh data
   * is exchanged during the @p initialize() call.
   *
   * @see getMeshVertexSize() to get the amount of vertices in the mesh
   * @see getMeshDimensions() to get the spacial dimensionality of the mesh
   */
  void getMeshVertexIDsAndCoordinates(
      ::precice::string_view    meshName,
      ::precice::span<VertexID> ids,
      ::precice::span<double>   coordinates) const;

  ///@}

  /** @name Experimental: Gradient Data
   * These API functions are \b experimental and may change in future versions.
   */
  ///@{

  /**
   * @brief Checks if the given data set requires gradient data.
   * We check if the data object has been initialized with the gradient flag.
   *
   * @experimental
   *
   * preCICE may require gradient data information from the solver and
   * ignores any API calls regarding gradient data if it is not required.
   * (When applying a nearest-neighbor-gradient mapping)
   *
   * @param[in] meshName the name of mesh that hold the data.
   * @param[in] dataName the name of the data.
   * @returns whether gradient is required
   */
  bool requiresGradientDataFor(::precice::string_view meshName,
                               ::precice::string_view dataName) const;

  /**
   * @brief Writes vector gradient data to a mesh.
   *
   * @experimental
   *
   * This function writes gradient values of specified vertices to data of a mesh.
   * Values are provided as a block of continuous memory defined by gradients.
   * The order of the provided gradient data follows the order specified by vertices.
   *
   * Each gradient or Jacobian depends on the dimensionality of the mesh and data.
   * Each gradient has a total of \ref getMeshDimensions() "getMeshDimensions(meshName)" * \ref getDataDimensions() "getDataDimensions(meshName, dataName)" components and
   * is stored in a linearised format as follows:
   *
   * | Spatial Dimensions | Scalar Data | Vectorial Data |
   * | --- | --- | --- |
   * | **2D** | s dx, s dy | x dx, y dx, x dy, y dy |
   * | **3D** | s dy, s dy, s dz | x dx, y dx, z dx, x dy, y dy, z dy, x dz, y dz, z dz |
   *
   *
   * The gradients/Jacobian for all \ref vertices are then contiguously saved in memory.
   *
   * Example for 2D Vectorial:
   *
   * | Index     | 0 | 1 | 2 | 3 | ... | 4n | 4n+1 | 4n+2 | 4n+3 |
   * | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
   * | Component | x0 dx| y0 dx | x0 dy | y0 dy | ... | xn dx | yn dx | xn dy | yn dy |
   *
   * @param[in] meshName the name of mesh that hold the data.
   * @param[in] dataName the name of the data to write to.
   * @param[in] vertices the vertex ids of the vertices to write gradient data to.
   * @param[in] gradients the linearised gradient data to write to preCICE.
   *
   * @pre Data has attribute hasGradient = true
   * @pre every VertexID in vertices is a return value of setMeshVertex or setMeshVertices
   * @pre gradients.size() == vertices.size() * getMeshDimensions(meshName) * getDataDimensions(meshName, dataName)
   *
   * @see Participant::setMeshVertex()
   * @see Participant::setMeshVertices()
   * @see Participant::getMeshDimensions()
   * @see Participant::getDataDimensions()
   *
   */
  void writeGradientData(
      ::precice::string_view          meshName,
      ::precice::string_view          dataName,
      ::precice::span<const VertexID> vertices,
      ::precice::span<const double>   gradients);

  ///@}

  /// Disable copy construction
  Participant(const Participant &copy) = delete;

  /// Disable assignment construction
  Participant &operator=(const Participant &assign) = delete;

private:
  /// Pointer to implementation of Participant.
  std::unique_ptr<impl::ParticipantImpl> _impl;

  // @brief To allow white box tests.
  friend struct testing::WhiteboxAccessor;
};

} // namespace precice
