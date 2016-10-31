/*!
* \file precice.hpp
* \brief Header of adapter class for coupling SU2 with preCICE for FSI.
* \author Alexander Rusch
*/



#include "precice/SolverInterface.hpp"
#include "SU2_CFD.hpp"
#include <string>
#include <stdlib.h>
using namespace precice;
using namespace std;

/*!
* \class Precice
* \brief Main class for coupling SU2 with preCICE.
*/
class Precice {
private:
  //MPI-rank(=index) and size
  int solverProcessIndex, solverProcessSize;
  //Coupling interface object
  SolverInterface solverInterface;
  //Mesh and boundary information
  CGeometry*** geometry_container;
  //Solution information
  CSolver**** solver_container;
  //Configuration information
  CConfig** config_container;
  //Grid movement information
  CVolumetricMovement** grid_movement;

  //General preCICE-related Variables
  int nDim;
  unsigned long *vertexSize; //Number of nodes at the wet surface (for each wet surface)
  short *valueMarkerWet; //List of wet surface marker values
  int **vertexIDs;
  int *displID;
  int *forceID;
  int *displDeltaID;
  double *forces;
  double *displacements;
  double *displacementDeltas;
  const string& coric;
  const string& cowic;
  bool processWorkingOnWetSurface;
  bool verbosityLevel_high;
  unsigned long globalNumberWetSurfaces; //Number of wet surfaces of the whole problem
  unsigned long localNumberWetSurfaces; //Number of wet surfaces, which this process is actually working on
  short *indexMarkerWetMappingLocalToGlobal; //Mapping relations for wet surfaces

  //Variables for storing the old state to reset to in case of an implicit simulation
  int nPoint; //Overall number of nodes of the problem
  int nVar; //Number of variables of the problem
  double **Coord_Saved, **Coord_n_Saved, **Coord_n1_Saved, **Coord_p1_Saved, **GridVel_Saved, ***GridVel_Grad_Saved;
  double dt_savedState;
  bool StopCalc_savedState;
  double **solution_Saved, **solution_time_n_Saved, **solution_time_n1_Saved;


public:
  /*!
  * \brief Constructor of the class.
  * \param[in] solverProcessIndex - Index of this solver Process.
  * \param[in] solverProcessSize - Overall number of solver Processes.
  * \param[in] geometry_container - Contains all geometric information of SU2.
  * \param[in] solver_container - Contains all solution and solver information of SU2.
  * \param[in] config_container - Contains all configuration information of SU2.
  * \param[in] grid_movement - Contains all grid movement information of SU2.
  */
  Precice( int solverProcessIndex, int solverProcessSize, CGeometry*** geometry_container, CSolver**** solver_container, CConfig** config_container, CVolumetricMovement** grid_movement );

  /*!
  * \brief Destructor of the class.
  */
  ~Precice();

  /*!
  * \brief Configures preCICE from the given xml file.
  *
  * Only after the configuration a usable state of a SolverInterface
  * object is achieved. However, most of the functionalities in preCICE can be
  * used only after initialize() has been called. Some actions, e.g. specifying
  * the solvers interface mesh, have to be done before initialize is called.
  *
  * In configure, the following is done:
  * - The XML configuration for preCICE is parsed and all objects containing
  *   data are created, but not necessarily filled with data.
  * - If a server is used, communication to that server is established.
  *
  * \param[in] configurationFileName - Name (with path) of the xml configuration
  *        file to be read.
  */
  void configure( const string& preciceConfigurationFileName );

  /*!
  * \brief Fully initializes preCICE to be used.
  *
  * Preconditions:
  * - configure() has been called successfully.
  *
  * Postconditions:
  * - Communication to the coupling partner/s is setup.
  * - Meshes are created from geometries or sent/received to/from coupling partners.
  * - If the solver is not starting the simulation, coupling data is received
  *   from the coupling partner's first computation.
  * - The length limitation of the first solver timestep is computed and returned.
  *
  * \return Maximum length of first timestep to be computed by the solver.
  */
  double initialize();

  /*!
  * \brief Advances preCICE after the solver has computed one timestep.
  *
  * Preconditions:
  * - initialize() has been called successfully.
  * - The solver calling advance() has computed one timestep.
  * - The solver has read and written all coupling data.
  * - The solver has the length of the timestep used available to be passed to
  *   preCICE.
  * - finalize() has not yet been called.
  *
  * Postconditions:
  * - Coupling data values specified to be exchanged in the configuration are
  *   exchanged with coupling partner/s.
  * - Coupling scheme state (computed time, computed timesteps, ...) is updated.
  * - For staggered coupling schemes, the coupling partner has computed one
  *   timestep/iteration with the coupling data written by this participant
  *   available.
  * - The coupling state is printed.
  * - Meshes with data are exported to files if specified in the configuration.
  * - The length limitation of the next solver timestep is computed and returned.
  *
  * \param[in] computedTimestepLength - Length of timestep computed by solver.
  * \return Maximum length of next timestep to be computed by solver.
  */
  double advance( double computedTimestepLength );

  /*!
  * \brief Checks wether preCICE coupling is still ongoing or not.
  *
  * \return True if coupling is still ongoing or false if coupling has ended.
  */
  bool isCouplingOngoing();

  /*!
  * \brief Checks wether an action (read or write iteration checkpoint) is required or not.
  *
  * \param[in] action - Reference to the strings coric or cowic, which determine wether read (coric) or write (cowic) action is checked.
  * \return True if an action is required or false if no action is required.
  */
  bool isActionRequired( const string& action );

  /*!
  * \brief Get reference to the string cowic which determines wether an iteration checkpoint needs to be written or not.
  *
  * \return The aforementioned reference to the cowic string.
  */
  const string& getCowic();

  /*!
  * \brief Get reference to the string coric which determines wether an iteration checkpoint needs to be read or not.
  *
  * \return The aforementioned reference to the coric string.
  */
  const string& getCoric();

  /*!
  * \brief Save the current state before the next SU2 run
  *
  * \param[in] StopCalc - Pointer to the bool StopCalc, which determines wether SU2 finishes computation after the current run or not. This variable must be passed so that its value can be stored and reloaded later.
  * \param[in] double - Pointer to the double dt, which is the current SU2 time step size. This variable must be passed so that its value can be stored and reloaded later.
  */
  void saveOldState( bool *StopCalc, double *dt );

  /*!
  * \brief Reload the state before the current SU2 run
  *
  * \param[in] StopCalc - Pointer to the bool StopCalc, which determines wether SU2 finishes computation after the current run or not. This variable must be passed so that its value can be reloaded from the saved state.
  * \param[in] double - Pointer to the double dt, which is the current time step size. This variable must be passed so that its value can be reloaded from the saved state.
  */
  void reloadOldState( bool *StopCalc, double *dt );

  /*!
  * \brief Finalizes preCICE.
  *
  * Preconditions:
  * - initialize() has been called successfully.
  * - isCouplingOngoing() has returned false.
  *
  * Postconditions:
  * - Communication channels are closed.
  * - Meshes and data are deallocated.
  */
  void finalize();

};
