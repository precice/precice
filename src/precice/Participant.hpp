#pragma once

#include <memory>
#include <precice/Version.h>
#include <precice/export.h>
#include <precice/span.hpp>
#include "precice/Types.hpp"

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

/** @brief forwards-compatible typedef for C++17 \ref std::string_view.
 *
 * The \ref string_view is a read-only view into sequence of characters.
 *
 * A good mental model for this type is a struct containing a pointer and a size like the following:
 *
 * @code{cpp}
 * struct string_view {
 *   const char* first;
 *   size_t size;
 * };
 * @endcode
 *
 * It can be constructed from:
 * - a pointer to a NULL-terminated C-string
 * - a pointer and a size
 * - a pointer to the first and one past the last element
 * - \ref std::string can directly convert itself to \ref std::string_view.
 *
 * In practise, using a function such as \ref Participant::getMeshDimensions() that expects a string_view looks like this:
 *
 * You can directly pass the string literal:
 *
 * @code{.cpp}
 * getMeshDimensions("MyMesh");
 * @endcode
 *
 * Or you can pass something string-like:
 *
 * @code{.cpp}
 * // one of
 * const char*      meshName = "MyMesh";
 * std::string      meshName = "MyMesh";
 * // followed by
 * getMeshDimensions(meshName);
 * @endcode
 *
 */
using string_view = ::precice::span<const char>;

/**
 * @brief Main Application Programming Interface of preCICE. Include using  `#include <precice/precice.hpp>`.
 *
 * @hidecollaborationgraph
 *
 * The flow of the API looks as follows:
 *
 * @startuml
 * skinparam conditionEndStyle hline
 * start
 * :Participant();
 * note right: Create a participant
 * :setMeshVertices();
 * note right: Define your meshes
 * rectangle {
 * note right
 * //Define mesh connectivity//
 * ----
 * Only required for
 * * projection mappings
 * * cell mappings
 * * watch integrals
 * * scaled-consistent mappings
 * end note
 * if ( requiresMeshConnectivityFor()) then (yes)
 * :setMeshEdges()
 * setMeshTriangles()
 * setMeshTetrahedra();
 * else (no)
 * endif
 * }
 *
 * rectangle {
 *
 * note right
 * //Provide initial data//
 * ----
 * Only required for non-zero
 * boundary conditions.
 * end note
 *
 * if (requiresInitialData()) then (yes)
 * :writeData();
 * else (no)
 * endif
 * }
 *
 * :initialize();
 * note right: Initialize coupling
 *
 * while (isCouplingOnGoing()) is (yes)
 *
 * rectangle {
 * note right
 * //Implicit coupling//
 * ----
 * New time window
 * Save solver state
 * end note
 * if (requiresWritingCheckpoint()) then (yes)
 * :solver writes checkpoint;
 * else (no)
 * endif
 * }
 *
 * :precice_dt = getMaxTimeStepSize()
 * solver_dt = solverGetAdaptiveDt()
 * dt = min(precice_dt, solver_dt);
 * note right: Agree on time step size
 *
 * :readData()
 * solverDoTimeStep(dt)
 * writeData()
 * advance(dt);
 * note right: Compute time step
 *
 * rectangle {
 * note right
 * //Implicit coupling//
 * ----
 * Iteration didn't converge
 * Restore solver state
 * ----
 * Iteration converged
 * Move solver to next time window
 * end note
 * if (requiresReadingCheckpoint()) then (yes)
 * :solver reads checkpoint;
 * else (no)
 * :solver moves in time;
 * endif
 * }
 *
 * endwhile (no)
 *
 * stop
 * @enduml
 *
 *
 *  @note
 *  We use solver, simulation code, and participant as synonyms.
 *  The preferred name in the documentation is participant.
 */
class PRECICE_API Participant {
public:
  /**
   * @name Construction and Configuration
   *
   * The API of preCICE is accessible via the \ref Participant class.
   * A constructed \ref Participant directly corresponds to a participant in the configuration file.
   *
   * Constructors require defining the parallel context of the Participant by providing
   * index and size of the current process, which are equivalent to rank and size in the MPI terminology.
   *
   * If preCICE is compiled with MPI, then there are multiple ways for it to be used:
   * 1. if a custom communicator is provided, preCICE uses it
   * 2. if MPI is already initialized, preCICE uses the MPI_COMM_WORLD
   * 3. otherwise, preCICE initializes MPI itself and uses the MPI_COMM_WORLD
   *
   * The MPI initialization is independent of which kind IntraCommunicator is configured.
   *
   * @{
   */

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

  /** @name Steering Methods
   *
   * The steering methods are responsible for the coupling logic in preCICE.
   *
   * 1. \ref initialize() connects participants, handles partitioned meshes, initializes coupling data, and computes mappings.
   * 2. \ref advance() steps the participant forwards in time, exchanging data, applying acceleration,
   * and transparently handling subcycling as well as implicit coupling.
   * 3. \ref finalize() closes down communication and optionally waits for other participants.
   * This is implicitly called in \ref ~Participant().
   *
   * @{
   */

  /**
   * @brief Fully initializes preCICE and coupling data.
   *
   * - Sets up a connection to the other participants of the coupled simulation.
   * - Pre-processes defined meshes and handles partitions in parallel.
   * - Receives first coupling data. The starting values for coupling data are zero by default.
   * - Determines maximum allowed size of the first time step to be computed.
   *
   * @pre initialize() has not yet been called.
   *
   * @post Parallel communication to the coupling partner(s) is setup.
   * @post Meshes are exchanged between coupling partners and the parallel partitions are created.
   * @post Initial coupling data was exchanged.
   *
   * @see getMaxTimeStepSize()
   * @see requiresInitialData()
   */
  void initialize();

  /**
   * @brief Advances preCICE after the solver has computed one time step.
   *
   * There are two cases to distinguish at this point.
   * If \p computedTimeStepSize == \ref getMaxTimeStepSize(), then the solver has reached
   * the end of a time window and proceeds the coupling.
   * This call is a computationally expensive process as it involves among other tasks:
   *
   * - Sending and resetting coupling data written by solver to coupling partners.
   * - Receiving coupling data read by solver.
   * - Computing and applying data mappings.
   * - Computing convergence measures in implicit coupling schemes.
   * - Computing acceleration of coupling data.
   * - Exchanging and computing information regarding the state of the coupled simulation.
   *
   * If \p computedTimeStepSize < \ref getMaxTimeStepSize(), then the solver hasn't reached the end of a time window and it is subcycling.
   * Depending on the configuration, written data can be used by preCICE to generate additional samples allowing for time interpolation using \ref readData().
   * This call is computationally inexpensive.
   *
   * @param[in] computedTimeStepSize Size of time step used by the solver.
   *
   * @attention All ranks of participants running in parallel \b must pass the same \p computedTimeStep to advance().
   *
   * @see getMaxTimeStepSize to get the maximum allowed value for \p computedTimeStepSize.
   *
   * @pre initialize() has been called successfully.
   * @pre The solver has computed one time step.
   * @pre The solver has written all coupling data.
   * @pre isCouplingOngoing() returns true.
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

  /** @name Implicit Coupling
   *
   * These functions are only required when you configure implicit coupling schemes.
   *
   * Implicitly-coupled solvers may iterate each time step until convergence is achieved.
   * This generally results in an outer loop handling time steps and an inner loop handling iterations.
   *
   * The preCICE abstracts the inner loop away by managing the iteration logic in advance.
   * If implicit coupling is configured, the adapter is required to either read or
   * write checkpoints using the API \ref requiresWritingCheckpoint() and \ref requiresReadingCheckpoint().
   *
   * The general flow looks as follows:
   *
   * @startuml
   * skinparam ConditionEndStyle hline
   * start
   * :initialize();
   *
   * while (isCouplingOngoing()) is (yes)
   *
   * if (requiresWritingCheckpoint()) then (yes)
   * :save solver checkpoint;
   * else (no)
   * endif
   *
   * :readData()
   * solve time step
   * writeData()
   * advance();
   *
   * if (requiresReadingCheckpoint()) then (yes)
   * : restore solver checkpoint;
   * else (no)
   * :move to next time window;
   * endif
   *
   *
   * endwhile (no)
   *
   * stop
   *
   * @enduml
   *
   *
   * @{
   */

  /** Checks if the participant is required to write an iteration checkpoint.
   *
   * If true, the participant is required to write an iteration checkpoint before calling advance().
   *
   * @note If implicit coupling is configured for this Participant, then this function **needs** to be called.
   *
   * @pre initialize() has been called
   *
   * @see requiresReadingCheckpoint()
   */
  bool requiresWritingCheckpoint();

  /** Checks if the participant is required to read an iteration checkpoint.
   *
   * If true, the participant is required to read an iteration checkpoint before calling advance().
   *
   * @note If implicit coupling is configured for this Participant, then this function **needs** to be called.
   *
   * @note This function returns false before the first call to advance().
   *
   * @pre initialize() has been called
   *
   * @see requiresWritingCheckpoint()
   */
  bool requiresReadingCheckpoint();

  ///@}

  /**
   * @name Status Queries
   *
   * In addition to the steering methods, the participant still needs some information regarding the coupling state.
   * Use \ref isCouplingOngoing() to check if the simulation has reached its end.
   * Control your time stepping using \ref getMaxTimeStepSize().
   *
   * To correctly size input and output buffers, you can access the dimensionality
   * of meshes using \ref getMeshDimensions() and data using \ref getDataDimensions().
   *
   * @{
   */

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
   * @brief Removes all vertices and connectivity information from the mesh
   *
   * @experimental
   *
   * Allows redefining a mesh during runtime.
   * After the call to resetMesh(), the mesh vertices need to be set with setMeshVertex() and setMeshVertices() again.
   * Connectivity information may be set as well.
   *
   * Reading data from this mesh using readData() is not possible until the next call to advance().
   *
   * @param[in] meshName the name of the mesh to reset
   *
   * @pre initialize() has been called
   * @pre isCouplingOngoing() is true
   *
   * @post previously returned vertex ids from setMeshVertex() and setMeshVertices() of the given mesh are invalid.
   */
  void resetMesh(::precice::string_view meshName);

  /**
   * @brief Creates a mesh vertex
   *
   * @param[in] meshName the name of the mesh to add the vertex to.
   * @param[in] position the coordinates of the vertex.
   * @returns the id of the created vertex
   *
   * @pre either initialize() has not yet been called or resetMesh(meshName) has been called since the last call to initialize() or advance()
   * @pre position.size() == getMeshDimensions(meshName)
   *
   * @see getMeshDimensions()
   */
  VertexID setMeshVertex(
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
   * @param[in] coordinates a span to the coordinates of the vertices
   *            The 2D-format is (d0x, d0y, d1x, d1y, ..., dnx, dny)
   *            The 3D-format is (d0x, d0y, d0z, d1x, d1y, d1z, ..., dnx, dny, dnz)
   *
   * @param[out] ids The ids of the created vertices
   *
   * @pre either initialize() has not yet been called or resetMesh(meshName) has been called since the last call to initialize() or advance()
   * @pre \p coordinates.size() == getMeshDimensions(meshName) * ids.size()
   *
   * @see getDimensions()
   */
  void setMeshVertices(
      ::precice::string_view        meshName,
      ::precice::span<const double> coordinates,
      ::precice::span<VertexID>     ids);

  /**
   * @brief Sets a mesh edge from vertex IDs
   *
   * Ignored if preCICE doesn't require connectivity for the mesh.
   *
   * @note The order of vertices does not matter.
   *
   * @param[in] meshName name of the mesh to add the edge to
   * @param[in] first ID of the first vertex of the edge
   * @param[in] second ID of the second vertex of the edge
   *
   * @pre vertices with IDs first and second were added to the mesh with the name meshName
   */
  void setMeshEdge(
      ::precice::string_view meshName,
      VertexID               first,
      VertexID               second);

  /**
   * @brief Sets multiple mesh edges from vertex IDs
   *
   * vertices contain pairs of vertex indices for each edge to define.
   * The format follows: e1a, e1b, e2a, e2b, ...
   * Ignored if preCICE doesn't require connectivity for the mesh.
   *
   * @note The order of vertices per edge does not matter.
   *
   * @param[in] meshName the name of the mesh to add the n edges to
   * @param[in] ids an array containing 2n vertex IDs for n edges
   *
   * @pre vertices in \p ids were added to the mesh with the name meshName
   * @pre \p ids.size() is multiple of 2
   *
   * @see requiresMeshConnectivityFor()
   */
  void setMeshEdges(
      ::precice::string_view          meshName,
      ::precice::span<const VertexID> ids);

  /**
   * @brief Sets mesh triangle from vertex IDs.
   *
   *
   *
   * @note The order of vertices does not matter.
   *
   * @param[in] meshName name of the mesh to add the triangle to
   * @param[in] first ID of the first vertex of the triangle
   * @param[in] second ID of the second vertex of the triangle
   * @param[in] third ID of the third vertex of the triangle
   *
   * @pre vertices with IDs first, second, and third were added to the mesh with the name meshName
   *
   * @see requiresMeshConnectivityFor()
   */
  void setMeshTriangle(
      ::precice::string_view meshName,
      VertexID               first,
      VertexID               second,
      VertexID               third);

  /**
   * @brief Sets multiple mesh triangles from vertex IDs
   *
   * vertices contain triples of vertex indices for each triangle to define.
   * The format follows: t1a, t1b, t1c, t2a, t2b, t2c, ...
   * Ignored if preCICE doesn't require connectivity for the mesh.
   *
   * @note The order of vertices per triangle does not matter.
   *
   * @param[in] meshName name of the mesh to add the n triangles to
   * @param[in] ids an array containing 3n vertex IDs for n triangles
   *
   * @pre vertices in \p ids were added to the mesh with the name meshName
   * @pre \p ids.size() is multiple of 3
   *
   * @see requiresMeshConnectivityFor()
   */
  void setMeshTriangles(
      ::precice::string_view          meshName,
      ::precice::span<const VertexID> ids);

  /**
   * @brief Sets a planar surface mesh quadrangle from vertex IDs.
   *
   * The planar quad will be triangulated, maximizing area-to-circumference.
   * Ignored if preCICE doesn't require connectivity for the mesh.
   *
   * @warning The order of vertices does not matter, however, only planar quads are allowed.
   *
   * @param[in] meshName name of the mesh to add the Quad to
   * @param[in] first ID of the first vertex of the Quad
   * @param[in] second ID of the second vertex of the Quad
   * @param[in] third ID of the third vertex of the Quad
   * @param[in] fourth ID of the fourth vertex of the Quad
   *
   * @pre vertices with IDs first, second, third, and fourth were added to the mesh with the name meshName
   *
   * @see requiresMeshConnectivityFor()
   */
  void setMeshQuad(
      ::precice::string_view meshName,
      VertexID               first,
      VertexID               second,
      VertexID               third,
      VertexID               fourth);

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
   * @param[in] meshName name of the mesh to add the n quads to
   * @param[in] ids an array containing 4n vertex IDs for n quads
   *
   * @pre vertices in \p ids were added to the mesh with the name meshName
   * @pre \p ids.size() is multiple of 4
   *
   * @see requiresMeshConnectivityFor()
   */
  void setMeshQuads(
      ::precice::string_view          meshName,
      ::precice::span<const VertexID> ids);

  /**
   * @brief Set tetrahedron in 3D mesh from vertex ID
   *
   * @note The order of vertices does not matter.
   *
   * @param[in] meshName name of the mesh to add the Tetrahedron to
   * @param[in] first ID of the first vertex of the Tetrahedron
   * @param[in] second ID of the second vertex of the Tetrahedron
   * @param[in] third ID of the third vertex of the Tetrahedron
   * @param[in] fourth ID of the fourth vertex of the Tetrahedron
   *
   * @pre vertices with IDs first, second, third, and fourth were added to the mesh with the name meshName
   *
   * @see requiresMeshConnectivityFor()
   */
  void setMeshTetrahedron(
      ::precice::string_view meshName,
      VertexID               first,
      VertexID               second,
      VertexID               third,
      VertexID               fourth);

  /**
   * @brief Sets multiple mesh tetrahedra from vertex IDs
   *
   * vertices contain quadruples of vertex indices for each tetrahedron to define.
   * The format follows: t1a, t1b, t1c, t1d, t2a, t2b, t2c, t2d, ...
   * Ignored if preCICE doesn't require connectivity for the mesh.
   *
   * @note The order of vertices per tetrahedron does not matter.
   *
   * @param[in] meshName name of the mesh to add the n tetrahedra to
   * @param[in] ids an array containing 4n vertex IDs for n tetrahedra
   *
   * @pre vertices in \p ids were added to the mesh with the name meshName
   * @pre ids.size() is multiple of 4
   *
   * @see requiresMeshConnectivityFor()
   */
  void setMeshTetrahedra(
      ::precice::string_view          meshName,
      ::precice::span<const VertexID> ids);

  ///@}

  /**
   * @name Data Access
   *
   * Data in preCICE is always associated to vertices on a \ref precice-mesh-access "defined mesh".
   * Use \ref getDataDimensions() to get the dimensionality of a data on a mesh.
   *
   * In each time step, you can access data on a mesh using \ref writeData() and \ref readData().
   * Calling \ref advance() may use written data to create a new sample in time, maps data between meshes, and communicates between participants.
   *
   * If you perform multiple time steps per time window, then preCICE may decide to keep samples of written data
   * to enable configured higher-order time interpolation in coupled participants.
   * The time interpolation is implemented by the relative time in \ref readData().
   * Written data is reset to 0 after each call to \ref advance().
   *
   * All data is initialized to 0 by default.
   * If you configure preCICE to provide custom initial data, then participants need to provide this data before calling \ref initialize().
   * After you defined the meshes, use \ref requiresInitialData() to check if initial data is required.
   * Then use \ref writeData() to specify your initial data and continue to \ref initialize().
   *
   * @{
   */

  /** Checks if the participant is required to provide initial data.
   *
   * If true, then the participant needs to write initial data to defined vertices
   * prior to calling initialize().
   *
   * @note If initial data is configured, then this function **needs** to be called.
   *
   * @pre initialize() has not yet been called
   */
  bool requiresInitialData();

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
   * @param[in] ids the vertex ids of the vertices to write data to.
   * @param[in] values the values to write to preCICE.
   *
   * @pre every VertexID in \p ids is a return value of setMeshVertex or setMeshVertices
   * @pre values.size() == getDataDimensions(meshName, dataName) * ids.size()
   *
   * @see Participant::setMeshVertex()
   * @see Participant::setMeshVertices()
   * @see Participant::getDataDimensions()
   */
  void writeData(
      ::precice::string_view          meshName,
      ::precice::string_view          dataName,
      ::precice::span<const VertexID> ids,
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
   * @param[in] ids the vertex ids of the vertices to read data from.
   * @param[in] relativeReadTime Point in time where data is read relative to the beginning of the current time step.
   * @param[out] values the destination memory to read the data from.
   *
   * @pre every VertexID in ids is a return value of setMeshVertex or setMeshVertices
   * @pre values.size() == getDataDimensions(meshName, dataName) * ids.size()
   * @pre resetMesh(meshName) has not been called since the last call to Participant::initialize() or Participant::advance()
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
      ::precice::span<const VertexID> ids,
      double                          relativeReadTime,
      ::precice::span<double>         values) const;
  ///@}

  /** @name Just-in-time mapping (or indirect access) (experimental)
   *
   * If one of your coupling meshes is not static and has varying locations over time, we can compute a data mapping
   * just-in-time. In such a case the user provides the coordinates of the moving mesh along with the API functions
   * \p mapAndWriteData or \p mapAndReadData to read and write data.
   *
   * The just-in-time mapping is closely connected to the \p Direct Access (see section below):
   *
   * Since one of the meshes is not given during the initialization, the user has to specify a region of interest
   * using \p setMeshAccessRegion() before calling \p initialize() to enable preCICE computing the repartitioning.
   *
   * Configuring this feature in the preCICE configuration file, requires two things:
   *
   * 1) The static mesh which is not moving (which is always a received mesh) needs api-access enabled
   * via `<receive-mesh name="StaticMesh" ... api-access="true"/>`. Similar to the Direct Access, the name of this static
   * mesh is then also the mesh name used in the API functions below, e.g., mapAndWriteData(StaticMesh, ...).
   *
   * 2) A mapping "from" or "to" the received mesh needs to be defined, where the "to" or "from" attribute in the configuration
   * needs to remain empty, e.g., ` <mapping:nearest-neighbor direction="read" from="StaticMesh" constraint="consistent" />`.
   * Note how the "to" attribute is not given in this configuration, as opposed to the conventional mapping configurations
   * in preCICE.
   *
   * @{
   */

  /**
   * @brief Writes data values to a mesh just-in-time (experimental).
   *
   * This function writes values of temporary vertices to data of a mesh.
   * As opposed to the writeData function using VertexIDs, this function allows to write data via coordinates,
   * which don't have to be specified during the initialization. This is particularly useful for meshes, which
   * vary over time. Note that using this function comes at a performance cost, since the specified mapping
   * needs to be computed locally for the given locations, whereas the other variant (writeData) can typically
   * exploit the static interface mesh and pre-compute data structures more efficient.
   *
   * Values are passed via a block of continuous memory defined by values in the order specified by vertices.
   *
   * The 1D/Scalar-format of values is (d0, d1, ..., dn)
   * The 2D-format of values is (d0x, d0y, d1x, d1y, ..., dnx, dny)
   * The 3D-format of values is (d0x, d0y, d0z, d1x, d1y, d1z, ..., dnx, dny, dnz)
   *
   * The time associated to the data is derived from the last call of \ref advance().
   *
   * @param[in] meshName the name of mesh that hold the data, typically a remote mesh from another participant.
   * @param[in] dataName the name of the data to write.
   * @param[in] coordinates a span to the coordinates of the vertices
   *        The 2D-format is (d0x, d0y, d1x, d1y, ..., dnx, dny)
   *        The 3D-format is (d0x, d0y, d0z, d1x, d1y, d1z, ..., dnx, dny, dnz)
   * @param[out] values the values containing the write data.
   *
   * @pre the coordinates are within the bounding box previously defined via \ref setMeshAccessRegion()
   *
   * @note the evaluated mapping computes the values corresponding to the initial configuration of the other provided mesh.
   * @note Only supported for conservative mapping constraints.
   * @note Caution when calling this function multiple times on the same data coordinates: There is no internal check and preCICE accumulates
   * data values for conservative mappings.
   * @note this function is currently part of the experimental API.
   *
   * @see Participant::setMeshAccessRegion()
   */
  void mapAndWriteData(
      ::precice::string_view        meshName,
      ::precice::string_view        dataName,
      ::precice::span<const double> coordinates,
      ::precice::span<const double> values);

  /**
   * @brief Reads data values from a mesh just-in-time. Values correspond to a given point in time relative to the beginning of the current timestep (experimental).
   *
   * This function reads values of temporary vertices from data of a mesh.
   * As opposed to the readData function using VertexIDs, this function allows to read data via coordinates,
   * which don't have to be specified during the initialization. This is particularly useful for meshes, which
   * vary over time. Note that using this function comes at a performance cost, since the specified mapping
   * needs to be computed locally for the given locations, whereas the other variant (readData) can typically
   * exploit the static interface mesh and pre-compute data structures more efficient.
   *
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
   * @param[in] meshName the name of mesh that hold the data, typically a remote mesh from another participant.
   * @param[in] dataName the name of the data to read from.
   * @param[in] coordinates a span to the coordinates of the vertices
   *        The 2D-format is (d0x, d0y, d1x, d1y, ..., dnx, dny)
   *        The 3D-format is (d0x, d0y, d0z, d1x, d1y, d1z, ..., dnx, dny, dnz)
   * @param[in] relativeReadTime Point in time where data is read relative to the beginning of the current time step.
   * @param[out] values the destination memory to read the data from.
   *
   * @pre the coordinates are within the bounding box previously defined via \ref setMeshAccessRegion()
   *
   * @post values contain the read data as specified in the above format.
   *
   * @note Note that the evaluated mapping computes the values corresponding to the initial configuration of the other provided mesh.
   *
   * @note this function is currently part of the experimental API.
   *
   * @see Participant::setMeshAccessRegion()
   */
  void mapAndReadData(
      ::precice::string_view        meshName,
      ::precice::string_view        dataName,
      ::precice::span<const double> coordinates,
      double                        relativeReadTime,
      ::precice::span<double>       values) const;

  ///@}

  /** @name Direct Access
   *
   * If you want or need to provide your own data mapping scheme, then you
   * can use direct mesh access to directly modify data on a received mesh.
   *
   * This requires to specify a region of interest using \ref setMeshAccessRegion() before calling \ref initialize().
   *
   * After \ref initialize(), you can use \ref getMeshVertexIDsAndCoordinates() to receive information on the received mesh.
   * Use the coordinates to compute your own data mapping scheme, and use the vertex IDs to read data form and write data to the mesh.
   *
   * @{
   */

  /**
   * @brief setMeshAccessRegion Define a region of interest on a received mesh
   *        (<receive-mesh ... from="otherParticipant />") in order to receive
   *        only a certain mesh region. Have a look at the website under
   *        https://precice.org/couple-your-code-direct-access.html or
   *        navigate manually to the page  Docs->Couple your code
   *        -> Advanced topics -> Accessing received meshes directly for
   *        a comprehensive documentation
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
   * @note This function can only be called once per mesh and rank
   * and trying to call it more than once results in an error.
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
   * @param[in] ids the vertex ids of the vertices to write gradient data to.
   * @param[in] gradients the linearised gradient data to write to preCICE.
   *
   * @pre Data has attribute hasGradient = true
   * @pre every VertexID in vertices is a return value of setMeshVertex or setMeshVertices
   * @pre gradients.size() == ids.size() * getMeshDimensions(meshName) * getDataDimensions(meshName, dataName)
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
      ::precice::span<const VertexID> ids,
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
