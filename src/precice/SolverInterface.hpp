#pragma once

#include <memory>
#include <set>
#include <string>
#include <vector>
#include "precice/Version.h"
#include "precice/export.h"

/**
 * forward declarations.
 */
namespace precice {
namespace impl {
class SolverInterfaceImpl;
}
namespace testing {
struct WhiteboxAccessor;
}
} // namespace precice

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {

/**
 * @brief Main Application Programming Interface of preCICE
 *
 * To adapt a solver to preCICE, follow the following main structure:
 *
 * -# Create an object of SolverInterface with SolverInterface()
 * -# Initialize preCICE with SolverInterface::initialize()
 * -# Advance to the next (time)step with SolverInterface::advance()
 * -# Finalize preCICE with SolverInterface::finalize()
 *
 *  @note
 *  We use solver, simulation code, and participant as synonyms.
 *  The preferred name in the documentation is participant.
 */
class PRECICE_API SolverInterface {
public:
  ///@name Construction and Configuration
  ///@{

  /**
   * @brief Constructs a SolverInterface for the given participant
   *
   * @param[in] participantName Name of the participant using the interface. Has to
   *        match the name given for a participant in the xml configuration file.
   * @param[in] configurationFileName Name (with path) of the xml configuration file.
   * @param[in] solverProcessIndex If the solver code runs with several processes,
   *        each process using preCICE has to specify its index, which has to start
   *        from 0 and end with solverProcessSize - 1.
   * @param[in] solverProcessSize The number of solver processes using preCICE.
   */
  SolverInterface(
      const std::string &participantName,
      const std::string &configurationFileName,
      int                solverProcessIndex,
      int                solverProcessSize);

  /**
   * @brief Constructs a SolverInterface for the given participant and a custom MPI communicator.
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
  SolverInterface(
      const std::string &participantName,
      const std::string &configurationFileName,
      int                solverProcessIndex,
      int                solverProcessSize,
      void *             communicator);

  ~SolverInterface();

  ///@}

  /// @name Steering Methods
  ///@{

  /**
   * @brief Fully initializes preCICE and coupling data.
   *
   * - Sets up a connection to the other participants of the coupled simulation.
   * - Creates all meshes, solver meshes need to be submitted before.
   * - Receives first coupling data. The starting values for coupling data are zero by default.
   * - Determines length of the first timestep to be computed.
   *
   * @pre initialize() has not yet been called.
   *
   * @post Parallel communication to the coupling partner(s) is setup.
   * @post Meshes are exchanged between coupling partners and the parallel partitions are created.
   * @post Initial coupling data was exchanged.
   *
   * @return Maximum length of first timestep to be computed by the solver.
   */
  double initialize();

  /**
   * @brief Advances preCICE after the solver has computed one timestep.
   *
   * - Sends and resets coupling data written by solver to coupling partners.
   * - Receives coupling data read by solver.
   * - Computes and applies data mappings.
   * - Computes acceleration of coupling data.
   * - Exchanges and computes information regarding the state of the coupled
   *   simulation.
   *
   * @param[in] computedTimestepLength Length of timestep used by the solver.
   *
   * @pre initialize() has been called successfully.
   * @pre The solver has computed one timestep.
   * @pre The solver has written all coupling data.
   * @pre isCouplngOngoing() returns true.
   * @pre finalize() has not yet been called.
   *
   * @post Coupling data values specified in the configuration are exchanged.
   * @post Coupling scheme state (computed time, computed timesteps, ...) is updated.
   * @post The coupling state is logged.
   * @post Configured data mapping schemes are applied.
   * @post [Second Participant] Configured acceleration schemes are applied.
   * @post Meshes with data are exported to files if configured.
   *
   * @return Maximum length of next timestep to be computed by solver.
   */
  double advance(double computedTimestepLength);

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
   * @brief Returns the number of spatial dimensions configured.
   *
   * @returns the configured dimension
   *
   * Currently, two and three dimensional problems can be solved using preCICE.
   * The dimension is specified in the XML configuration.
   */
  int getDimensions() const;

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
  //  void resetMesh ( int meshID );

  /**
   * @brief Checks if the mesh with given name is used by a solver.
   *
   * @param[in] meshName the name of the mesh
   * @returns whether the mesh is used.
   */
  bool hasMesh(const std::string &meshName) const;

  /**
   * @brief Checks if the given mesh requires connectivity.
   *
   * preCICE may require connectivity information from the solver and
   * ignores any API calls regarding connectivity if it is not required.
   * Use this function to conditionally generate this connectivity.
   *
   * @param[in] meshID the id of the mesh
   * @returns whether connectivity is required
   */
  bool requiresMeshConnectivityFor(int meshID) const;

  /**
   * @brief Creates a mesh vertex
   *
   * @param[in] meshID the id of the mesh to add the vertex to.
   * @param[in] position a pointer to the coordinates of the vertex.
   * @returns the id of the created vertex
   *
   * @pre initialize() has not yet been called
   * @pre count of available elements at position matches the configured dimension
   *
   * @see getDimensions()
   */
  int setMeshVertex(
      int           meshID,
      const double *position);

  /**
   * @brief Returns the number of vertices of a mesh.
   *
   * @param[in] meshID the id of the mesh
   * @returns the amount of the vertices of the mesh
   *
   * @pre This function can be called on received meshes as well as provided
   * meshes. However, you need to call this function after @p initialize(),
   * if the \p meshID corresponds to a received mesh, since the relevant mesh data
   * is exchanged during the @p initialize() call.
   */
  int getMeshVertexSize(int meshID) const;

  /**
   * @brief Creates multiple mesh vertices
   *
   * @param[in] meshID the id of the mesh to add the vertices to.
   * @param[in] size Number of vertices to create
   * @param[in] positions a pointer to the coordinates of the vertices
   *            The 2D-format is (d0x, d0y, d1x, d1y, ..., dnx, dny)
   *            The 3D-format is (d0x, d0y, d0z, d1x, d1y, d1z, ..., dnx, dny, dnz)
   *
   * @param[out] ids The ids of the created vertices
   *
   * @pre initialize() has not yet been called
   * @pre count of available elements at positions matches the configured dimension * size
   * @pre count of available elements at ids matches size
   *
   * @see getDimensions()
   */
  void setMeshVertices(
      int           meshID,
      int           size,
      const double *positions,
      int *         ids);

  /**
   * @brief Sets a mesh edge from vertex IDs
   *
   * Ignored if preCICE doesn't require connectivity for the mesh.
   *
   * @note The order of vertices does not matter.
   *
   * @param[in] meshID ID of the mesh to add the edge to
   * @param[in] firstVertexID ID of the first vertex of the edge
   * @param[in] secondVertexID ID of the second vertex of the edge
   *
   * @pre vertices with firstVertexID and secondVertexID were added to the mesh with the ID meshID
   */
  void setMeshEdge(
      int meshID,
      int firstVertexID,
      int secondVertexID);

  /**
   * @brief Sets multiple mesh edge from vertex IDs
   *
   * vertices contain pairs of vertex indices for each edge to define.
   * The format follows: e1a, e1b, e2a, e2b, ...
   * Ignored if preCICE doesn't require connectivity for the mesh.
   *
   * @note The order of vertices per edge does not matter.
   *
   * @param[in] meshID ID of the mesh to add the edges to
   * @param[in] size the amount of edges to set
   * @param[in] vertices an array containing 2*size vertex IDs
   *
   * @pre vertices were added to the mesh with the ID meshID
   *
   * @see requiresMeshConnectivityFor()
   */
  void setMeshEdges(
      int        meshID,
      int        size,
      const int *vertices);

  /**
   * @brief Sets mesh triangle from vertex IDs.
   *
   *
   *
   * @note The order of vertices does not matter.
   *
   * @param[in] meshID ID of the mesh to add the triangle to
   * @param[in] firstVertexID ID of the first vertex of the triangle
   * @param[in] secondVertexID ID of the second vertex of the triangle
   * @param[in] thirdVertexID ID of the third vertex of the triangle
   *
   * @pre edges with firstVertexID, secondVertexID, and thirdVertexID were added to the mesh with the ID meshID
   *
   * @see requiresMeshConnectivityFor()
   */
  void setMeshTriangle(
      int meshID,
      int firstVertexID,
      int secondVertexID,
      int thirdVertexID);

  /**
   * @brief Sets multiple mesh triangles from vertex IDs
   *
   * vertices contain triples of vertex indices for each triangle to define.
   * The format follows: t1a, t1b, t1c, t2a, t2b, t2c, ...
   * Ignored if preCICE doesn't require connectivity for the mesh.
   *
   * @note The order of vertices per triangle does not matter.
   *
   * @param[in] meshID ID of the mesh to add the triangles to
   * @param[in] size the amount of triangles to set
   * @param[in] vertices an array containing 3*size vertex IDs
   *
   * @pre vertices were added to the mesh with the ID meshID
   *
   * @see requiresMeshConnectivityFor()
   */
  void setMeshTriangles(
      int        meshID,
      int        size,
      const int *vertices);

  /**
   * @brief Sets a planar surface mesh quadrangle from vertex IDs.
   *
   * The planar quad will be triangulated, maximizing area-to-circumference.
   * Ignored if preCICE doesn't require connectivity for the mesh.
   *
   * @warning The order of vertices does not matter, however, only planar quads are allowed.
   *
   * @param[in] meshID ID of the mesh to add the Quad to
   * @param[in] firstVertexID ID of the first vertex of the Quad
   * @param[in] secondVertexID ID of the second vertex of the Quad
   * @param[in] thirdVertexID ID of the third vertex of the Quad
   * @param[in] fourthVertexID ID of the fourth vertex of the Quad
   *
   * @pre vertices with firstVertexID, secondVertexID, thirdVertexID, and fourthVertexID were added to the mesh with the ID meshID
   *
   * @see requiresMeshConnectivityFor()
   */
  void setMeshQuad(
      int meshID,
      int firstVertexID,
      int secondVertexID,
      int thirdVertexID,
      int fourthVertexID);

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
   * @param[in] meshID ID of the mesh to add the quad to
   * @param[in] size the amount of quads to set
   * @param[in] vertices an array containing 4*size vertex IDs
   *
   * @pre vertices were added to the mesh with the ID meshID
   *
   * @see requiresMeshConnectivityFor()
   */
  void setMeshQuads(
      int        meshID,
      int        size,
      const int *vertices);

  /**
   * @brief Set tetrahedron in 3D mesh from vertex ID
   *
   * @note The order of vertices does not matter.
   *
   * @param[in] meshID ID of the mesh to add the Tetrahedron to
   * @param[in] firstVertexID ID of the first vertex of the Tetrahedron
   * @param[in] secondVertexID ID of the second vertex of the Tetrahedron
   * @param[in] thirdVertexID ID of the third vertex of the Tetrahedron
   * @param[in] fourthVertexID ID of the fourth vertex of the Tetrahedron
   *
   * @pre vertices with firstVertexID, secondVertexID, thirdVertexID, and fourthVertexID were added to the mesh with the ID meshID
   *
   * @see requiresMeshConnectivityFor()
   */
  void setMeshTetrahedron(
      int meshID,
      int firstVertexID,
      int secondVertexID,
      int thirdVertexID,
      int fourthVertexID);

  /**
   * @brief Sets multiple mesh tetrahedra from vertex IDs
   *
   * vertices contain quadruples of vertex indices for each tetrahedron to define.
   * The format follows: t1a, t1b, t1c, t1d, t2a, t2b, t2c, t2d, ...
   * Ignored if preCICE doesn't require connectivity for the mesh.
   *
   * @note The order of vertices per tetrahedron does not matter.
   *
   * @param[in] meshID ID of the mesh to add the tetrahedra to
   * @param[in] size the amount of tetrahedra to set
   * @param[in] vertices an array containing 4*size vertex IDs
   *
   * @pre vertices were added to the mesh with the ID meshID
   *
   * @see requiresMeshConnectivityFor()
   */
  void setMeshTetrahedra(
      int        meshID,
      int        size,
      const int *vertices);

  ///@}

  ///@name Data Access
  ///@{

  /**
   * @brief Checks if the data with given name is used by a solver and mesh.
   *
   * @param[in] dataName the name of the data
   * @param[in] meshID the id of the associated mesh
   * @returns whether the mesh is used.
   */
  bool hasData(const std::string &dataName, int meshID) const;

  /**
   * @brief Writes vector data given as block.
   *
   * This function writes values of specified vertices to a dataID.
   * Values are provided as a block of continuous memory.
   * valueIndices contains the indices of the vertices
   *
   * The 2D-format of values is (d0x, d0y, d1x, d1y, ..., dnx, dny)
   * The 3D-format of values is (d0x, d0y, d0z, d1x, d1y, d1z, ..., dnx, dny, dnz)
   *
   * @param[in] dataID ID to write to.
   * @param[in] size Number n of vertices.
   * @param[in] valueIndices Indices of the vertices.
   * @param[in] values Pointer to the vector values.
   *
   * @pre count of available elements at values matches the configured dimension * size
   * @pre count of available elements at valueIndices matches the given size
   * @pre initialize() has been called
   *
   * @see SolverInterface::setMeshVertex()
   */
  void writeBlockVectorData(
      int           dataID,
      int           size,
      const int *   valueIndices,
      const double *values);

  /**
   * @brief Writes vector data to a vertex
   *
   * This function writes a value of a specified vertex to a dataID.
   * Values are provided as a block of continuous memory.
   *
   * The 2D-format of value is (x, y)
   * The 3D-format of value is (x, y, z)
   *
   * @param[in] dataID ID to write to.
   * @param[in] valueIndex Index of the vertex.
   * @param[in] value Pointer to the vector value.
   *
   * @pre count of available elements at value matches the configured dimension
   * @pre initialize() has been called
   *
   * @see SolverInterface::setMeshVertex()
   */
  void writeVectorData(
      int           dataID,
      int           valueIndex,
      const double *value);

  /**
   * @brief Writes scalar data given as block.
   *
   * This function writes values of specified vertices to a dataID.
   * Values are provided as a block of continuous memory.
   * valueIndices contains the indices of the vertices
   *
   * @param[in] dataID ID to write to.
   * @param[in] size Number n of vertices.
   * @param[in] valueIndices Indices of the vertices.
   * @param[in] values Pointer to the values.
   *
   * @pre count of available elements at values matches the given size
   * @pre count of available elements at valueIndices matches the given size
   * @pre initialize() has been called
   *
   * @see SolverInterface::setMeshVertex()
   */
  void writeBlockScalarData(
      int           dataID,
      int           size,
      const int *   valueIndices,
      const double *values);

  /**
   * @brief Writes scalar data to a vertex
   *
   * This function writes a value of a specified vertex to a dataID.
   *
   * @param[in] dataID ID to write to.
   * @param[in] valueIndex Index of the vertex.
   * @param[in] value The value to write.
   *
   * @pre initialize() has been called
   *
   * @see SolverInterface::setMeshVertex()
   */
  void writeScalarData(
      int    dataID,
      int    valueIndex,
      double value);

  /**
   * @brief Reads vector data values given as block from a mesh. Values correspond to the end of the current time window.
   *
   * This function reads values of specified vertices from a dataID.
   * Values are read into a block of continuous memory.
   * valueIndices contains the indices of the vertices.
   *
   * The 2D-format of values is (d0x, d0y, d1x, d1y, ..., dnx, dny)
   * The 3D-format of values is (d0x, d0y, d0z, d1x, d1y, d1z, ..., dnx, dny, dnz)
   *
   * @param[in] dataID ID to read from.
   * @param[in] size Number n of vertices.
   * @param[in] valueIndices Indices of the vertices.
   * @param[out] values Pointer to read destination.
   *
   * @pre count of available elements at values matches the configured dimension * size
   * @pre count of available elements at valueIndices matches the given size
   * @pre initialize() has been called
   *
   * @post values contain the read data as specified in the above format.
   *
   * @see SolverInterface::setMeshVertex()
   */
  void readBlockVectorData(
      int        dataID,
      int        size,
      const int *valueIndices,
      double *   values) const;

  /**
   * @brief Reads vector data at a vertex on a mesh. Values correspond to the end of the current time window.
   *
   * This function reads a value of a specified vertex from a dataID.
   * Values are provided as a block of continuous memory.
   *
   * The 2D-format of value is (x, y)
   * The 3D-format of value is (x, y, z)
   *
   * @param[in] dataID ID to read from.
   * @param[in] valueIndex Index of the vertex.
   * @param[out] value Pointer to the vector value.
   *
   * @pre count of available elements at value matches the configured dimension
   * @pre initialize() has been called
   *
   * @post value contains the read data as specified in the above format.
   *
   * @see SolverInterface::setMeshVertex()
   */
  void readVectorData(
      int     dataID,
      int     valueIndex,
      double *value) const;

  /**
   * @brief Reads scalar data values given as block from a mesh. Values correspond to the end of the current time window.
   *
   * This function reads values of specified vertices from a dataID.
   * Values are provided as a block of continuous memory.
   * valueIndices contains the indices of the vertices.
   *
   * @param[in] dataID ID to read from.
   * @param[in] size Number n of vertices.
   * @param[in] valueIndices Indices of the vertices.
   * @param[out] values Pointer to the read destination.
   *
   * @pre count of available elements at values matches the given size
   * @pre count of available elements at valueIndices matches the given size
   * @pre initialize() has been called
   *
   * @post values contains the read data.
   *
   * @see SolverInterface::setMeshVertex()
   */
  void readBlockScalarData(
      int        dataID,
      int        size,
      const int *valueIndices,
      double *   values) const;

  /**
   * @brief Reads scalar data at a vertex on a mesh. Values correspond to the end of the current time window.
   *
   * This function reads a value of a specified vertex from a dataID.
   *
   * @param[in] dataID ID to read from.
   * @param[in] valueIndex Index of the vertex.
   * @param[out] value Read destination of the value.
   *
   * @pre initialize() has been called
   *
   * @post value contains the read data.
   *
   * @see SolverInterface::setMeshVertex()
   */
  void readScalarData(
      int     dataID,
      int     valueIndex,
      double &value) const;

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
   * @param[in] meshID ID of the mesh you want to access through the bounding box
   * @param[in] boundingBox Axis aligned bounding boxes which has in 3D the format
   *            [x_min, x_max, y_min, y_max, z_min, z_max]
   *
   * @pre @p initialize() has not yet been called.
   */
  void setMeshAccessRegion(
      const int     meshID,
      const double *boundingBox) const;

  /**
   * @brief getMeshVerticesAndIDs Iterates over the region of
   *        interest defined by bounding boxes and reads the corresponding
   *        coordinates omitting the mapping.
   *
   * @experimental
   *
   * @param[in]  meshID corresponding mesh ID
   * @param[in]  size return value of @p getMeshVertexSize()
   * @param[out] ids ids corresponding to the coordinates
   * @param[out] coordinates the coordinates associated to the \p ids and
   *             corresponding data values (dim * \p size)
   *
   * @pre IDs and coordinates need to have the correct size, which can be queried by @p getMeshVertexSize()
   *
   * @pre This function can be called on received meshes as well as provided
   * meshes. However, you need to call this function after @p initialize(),
   * if the \p meshID corresponds to a received mesh, since the relevant mesh data
   * is exchanged during the @p initialize() call.
   */
  void getMeshVerticesAndIDs(
      const int meshID,
      const int size,
      int *     ids,
      double *  coordinates) const;

  ///@}

  /** @name Experimental: Time Interpolation
   * These API functions are \b experimental and may change in future versions.
   */
  ///@{

  /**
   * @brief Reads vector data values given as block from a mesh. Values correspond to a given point in time relative to the beginning of the current timestep.
   *
   * @experimental
   *
   * This function reads values of specified vertices from a dataID.
   * Values are read into a block of continuous memory.
   * valueIndices contains the indices of the vertices.
   *
   * The 2D-format of values is (d0x, d0y, d1x, d1y, ..., dnx, dny)
   * The 3D-format of values is (d0x, d0y, d0z, d1x, d1y, d1z, ..., dnx, dny, dnz)
   *
   * The data is read at relativeReadTime, which indicates the point in time measured from the beginning of the current time step.
   * relativeReadTime = 0 corresponds to data at the beginning of the time step. Assuming that the user will call advance(dt) at the
   * end of the time step, dt indicates the length of the current time step. Then relativeReadTime = dt corresponds to the data at
   * the end of the time step.
   *
   * @param[in] dataID ID to read from.
   * @param[in] size Number n of vertices.
   * @param[in] valueIndices Indices of the vertices.
   * @param[in] relativeReadTime Point in time where data is read relative to the beginning of the current time step.
   * @param[out] values Pointer to read destination.
   *
   * @pre count of available elements at values matches the configured dimension * size
   * @pre count of available elements at valueIndices matches the given size
   * @pre initialize() has been called
   *
   * @post values contain the read data as specified in the above format.
   *
   * @see SolverInterface::setMeshVertex()
   */
  void readBlockVectorData(
      int        dataID,
      int        size,
      const int *valueIndices,
      double     relativeReadTime,
      double *   values) const;

  /**
   * @brief Reads vector data at a vertex on a mesh. Values correspond to a given point in time relative to the beginning of the current timestep.
   *
   * @experimental
   *
   * This function reads a value of a specified vertex from a dataID.
   * Values are provided as a block of continuous memory.
   *
   * The 2D-format of value is (x, y)
   * The 3D-format of value is (x, y, z)
   *
   * The data is read at relativeReadTime, which indicates the point in time measured from the beginning of the current time step.
   * relativeReadTime = 0 corresponds to data at the beginning of the time step. Assuming that the user will call advance(dt) at the
   * end of the time step, dt indicates the length of the current time step. Then relativeReadTime = dt corresponds to the data at
   * the end of the time step.
   *
   * @param[in] dataID ID to read from.
   * @param[in] valueIndex Index of the vertex.
   * @param[in] relativeReadTime Point in time where data is read relative to the beginning of the current time step.
   * @param[out] value Pointer to the vector value.
   *
   * @pre count of available elements at value matches the configured dimension
   * @pre initialize() has been called
   *
   * @post value contains the read data as specified in the above format.
   *
   * @see SolverInterface::setMeshVertex()
   */
  void readVectorData(
      int     dataID,
      int     valueIndex,
      double  relativeReadTime,
      double *value) const;

  /**
   * @brief Reads scalar data values given as block from a mesh. Values correspond to a given point in time relative to the beginning of the current timestep.
   *
   * @experimental
   *
   * This function reads values of specified vertices from a dataID.
   * Values are provided as a block of continuous memory.
   * valueIndices contains the indices of the vertices.
   *
   * The data is read at relativeReadTime, which indicates the point in time measured from the beginning of the current time step.
   * relativeReadTime = 0 corresponds to data at the beginning of the time step. Assuming that the user will call advance(dt) at the
   * end of the time step, dt indicates the length of the current time step. Then relativeReadTime = dt corresponds to the data at
   * the end of the time step.
   *
   * @param[in] dataID ID to read from.
   * @param[in] size Number n of vertices.
   * @param[in] valueIndices Indices of the vertices.
   * @param[in] relativeReadTime Point in time where data is read relative to the beginning of the current time step.
   * @param[out] values Pointer to the read destination.
   *
   * @pre count of available elements at values matches the given size
   * @pre count of available elements at valueIndices matches the given size
   * @pre initialize() has been called
   *
   * @post values contains the read data.
   *
   * @see SolverInterface::setMeshVertex()
   */
  void readBlockScalarData(
      int        dataID,
      int        size,
      const int *valueIndices,
      double     relativeReadTime,
      double *   values) const;

  /**
   * @brief Reads scalar data at a vertex on a mesh. Values correspond to a given point in time relative to the beginning of the current timestep.
   *
   * @experimental
   *
   * This function reads a value of a specified vertex from a dataID.
   *
   * The data is read at relativeReadTime, which indicates the point in time measured from the beginning of the current time step.
   * relativeReadTime = 0 corresponds to data at the beginning of the time step. Assuming that the user will call advance(dt) at the
   * end of the time step, dt indicates the length of the current time step. Then relativeReadTime = dt corresponds to the data at
   * the end of the time step.
   *
   * @param[in] dataID ID to read from.
   * @param[in] valueIndex Index of the vertex.
   * @param[in] relativeReadTime Point in time where data is read relative to the beginning of the current time step
   * @param[out] value Read destination of the value.
   *
   * @pre initialize() has been called
   *
   * @post value contains the read data.
   *
   * @see SolverInterface::setMeshVertex()
   */
  void readScalarData(
      int     dataID,
      int     valueIndex,
      double  relativeReadTime,
      double &value) const;

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
   * @param[in] dataID the id of the data
   * @returns whether gradient is required
   */
  bool requiresGradientDataFor(int dataID) const;

  /**
   * @brief Writes vector gradient data given as block.
   *
   * @experimental
   *
   * This function writes values of specified vertices to a dataID.
   * Values are provided as a block of continuous memory.
   * \p valueIndices contains the indices of the vertices
   *
   * The values are passed in the same format applied in \ref writeVectorGradientData() for each data vertex:
   *
   * The 2D-format of \p gradientValues is ( v0x_dx, v0y_dx, v0x_dy, v0y_dy,
   *                                         v1x_dx, v1y_dx, v1x_dy, v1y_dy,
   *                                         ... ,
   *                                         vnx_dx, vny_dx, vnx_dy, vny_dy)
   *
   * corresponding to the vector data v0 = (v0x, v0y) , v1 = (v1x, v1y), ... , vn = (vnx, vny) differentiated in spatial directions x and y.
   *
   *
   * The 3D-format of \p gradientValues is ( v0x_dx, v0y_dx, v0z_dx, v0x_dy, v0y_dy, v0z_dy, v0x_dz, v0y_dz, v0z_dz,
   *                                         v1x_dx, v1y_dx, v1z_dx, v1x_dy, v1y_dy, v1z_dy, v1x_dz, v1y_dz, v1z_dz,
   *                                         ... ,
   *                                         vnx_dx, vny_dx, vnz_dx, vnx_dy, vny_dy, vnz_dy, vnx_dz, vny_dz, vnz_dz)
   *
   * corresponding to the vector data v0 = (v0x, v0y, v0z) , v1 = (v1x, v1y, v1z), ... , vn = (vnx, vny, vnz) differentiated in spatial directions x,y and z.
   *
   * @param[in] dataID ID to write to.
   * @param[in] size Number n of vertices.
   * @param[in] valueIndices Indices of the vertices.
   * @param[in] gradientValues Pointer to the gradient values.
   *
   * @pre count of available elements at gradient values matches the configured dimension * size
   * @pre count of available elements at valueIndices matches the given size
   * @pre initialize() has been called
   * @pre Data with dataID has attribute hasGradient = true
   *
   * @see SolverInterface::setMeshVertex()
   */
  void writeBlockVectorGradientData(
      int           dataID,
      int           size,
      const int *   valueIndices,
      const double *gradientValues);

  /**
   * @brief Writes scalar gradient data to a vertex
   *
   * @experimental
   *
   * This function writes a the corresponding gradient value of a specified vertex to a dataID.
   * Values are provided as a block of continuous memory.
   *
   * @param[in] dataID ID to write to.
   * @param[in] valueIndex Index of the vertex.
   * @param[in] gradientValues Gradient values differentiated in the spacial direction (dx, dy) for 2D space, (dx, dy, dz) for 3D space
   *
   * @pre count of available elements at value matches the configured dimension
   * @pre initialize() has been called
   * @pre vertex with dataID exists and contains data
   * @pre Data with dataID has attribute hasGradient = true
   *
   * @see SolverInterface::setMeshVertex()
   */
  void writeScalarGradientData(
      int           dataID,
      int           valueIndex,
      const double *gradientValues);

  /**
   * @brief Writes vector gradient data to a vertex
   *
   * @experimental
   *
   * This function writes the corresponding gradient matrix value of a specified vertex to a dataID.
   * Values are provided as a block of continuous memory.
   *
   * The gradients need to be provided in the following format:
   *
   * The 2D-format of \p gradientValues is (vx_dx, vy_dx, vx_dy, vy_dy) matrix corresponding to the data block v = (vx, vy)
   * differentiated respectively in x-direction dx and y-direction dy
   *
   * The 3D-format of \p gradientValues is (vx_dx, vy_dx, vz_dx, vx_dy, vy_dy, vz_dy, vx_dz, vy_dz, vz_dz) matrix
   * corresponding to the data block v = (vx, vy, vz) differentiated respectively in spatial directions x-direction dx and y-direction dy and z-direction dz
   *
   * @param[in] dataID ID to write to.
   * @param[in] valueIndex Index of the vertex.
   * @param[in] gradientValue pointer to the gradient value.
   *
   * @pre count of available elements at value matches the configured dimension
   * @pre initialize() has been called
   * @pre vertex with dataID exists and contains data
   * @pre Data with dataID has attribute hasGradient = true
   *
   * @see SolverInterface::setMeshVertex()
   */
  void writeVectorGradientData(
      int           dataID,
      int           valueIndex,
      const double *gradientValues);

  /**
   * @brief Writes scalar gradient data given as block.
   *
   * @experimental
   *
   * This function writes values of specified vertices to a dataID.
   * Values are provided as a block of continuous memory.
   * valueIndices contains the indices of the vertices
   *
   * The gradients need to be provided in the following format:
   *
   * The 2D-format of \p gradientValues is (v0_dx, v0_dy, v1_dx, v1_dy, ... , vn_dx, vn_dy, vn_dz)
   * corresponding to the scalar data v0, v1, ... , vn differentiated in spatial directions x and y.
   *
   * The 3D-format of \p gradientValues is (v0_dx, v0_dy, v0_dz, v1_dx, v1_dy, v1_dz, ... , vn_dx, vn_dy, vn_dz)
   * corresponding to the scalar data v0, v1, ... , vn differentiated in spatial directions x, y and z.
   *
   * @param[in] dataID ID to write to.
   * @param[in] size Number n of vertices.
   * @param[in] valueIndices Indices of the vertices.
   * @param[in] gradientValues Pointer to the gradient values.
   *
   * @pre count of available elements at values matches the given size
   * @pre count of available elements at valueIndices matches the given size
   * @pre initialize() has been called
   * @pre Data with dataID has attribute hasGradient = true
   *
   * @see SolverInterface::setMeshVertex()
   */
  void writeBlockScalarGradientData(
      int           dataID,
      int           size,
      const int *   valueIndices,
      const double *gradientValues);

  ///@}

  /// Disable copy construction
  SolverInterface(const SolverInterface &copy) = delete;

  /// Disable assignment construction
  SolverInterface &operator=(const SolverInterface &assign) = delete;

private:
  /// Pointer to implementation of SolverInterface.
  std::unique_ptr<impl::SolverInterfaceImpl> _impl;

  // @brief To allow white box tests.
  friend struct testing::WhiteboxAccessor;
};

} // namespace precice
