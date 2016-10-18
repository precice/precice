/*!
* \file precice.cpp
* \brief Adapter class for coupling SU2 with preCICE for FSI.
* \author Alexander Rusch
*/

#include "../include/precice.hpp"

Precice::Precice( int solverProcessIndex, int solverProcessSize, CGeometry*** geometry_container, CSolver**** solver_container, CConfig** config_container, CVolumetricMovement** grid_movement )
:
solverProcessIndex(solverProcessIndex),
solverProcessSize(solverProcessSize),
solverInterface( "SU2_CFD", solverProcessIndex, solverProcessSize ),
nDim(geometry_container[ZONE_0][MESH_0]->GetnDim()),
geometry_container(geometry_container),
solver_container(solver_container),
config_container(config_container),
grid_movement(grid_movement),
vertexIDs(NULL),
displID(NULL),
forceID(NULL),
displDeltaID(NULL),
forces(NULL),
displacementDeltas(NULL),
//For implicit coupling
coric(precice::constants::actionReadIterationCheckpoint()),
cowic(precice::constants::actionWriteIterationCheckpoint()),
processWorkingOnWetSurface(true),
verbosityLevel_high(config_container[ZONE_0]->GetpreCICE_VerbosityLevel_High()),
globalNumberWetSurfaces(config_container[ZONE_0]->GetpreCICE_NumberWetSurfaces()),
localNumberWetSurfaces(0),
//Get value (= index) of the marker corresponding to the wet surface
//It is implied, that only one marker is used for the entire wet surface, even if it is split into parts
valueMarkerWet(NULL),
vertexSize(NULL),
indexMarkerWetMappingLocalToGlobal(NULL),
//Variables for implicit coupling
nPoint(geometry_container[ZONE_0][MESH_0]->GetnPoint()),
nVar(solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetnVar()),
Coord_Saved(NULL),
Coord_n_Saved(NULL),
Coord_n1_Saved(NULL),
Coord_p1_Saved(NULL),
GridVel_Saved(NULL),
GridVel_Grad_Saved(NULL),
dt_savedState(0),
StopCalc_savedState(false),
solution_Saved(NULL),
solution_time_n_Saved(NULL),
solution_time_n1_Saved(NULL)
{
  Coord_Saved = new double*[nPoint];
  Coord_n_Saved = new double*[nPoint];
  Coord_n1_Saved = new double*[nPoint];
  Coord_p1_Saved = new double*[nPoint];
  GridVel_Saved = new double*[nPoint];
  GridVel_Grad_Saved = new double**[nPoint];
  solution_Saved = new double*[nPoint];
  solution_time_n_Saved = new double*[nPoint];
  solution_time_n1_Saved = new double*[nPoint];
  for (int iPoint = 0; iPoint < nPoint; iPoint++) {
    Coord_Saved[iPoint] = new double[nDim];
    Coord_n_Saved[iPoint] = new double[nDim];
    Coord_n1_Saved[iPoint] = new double[nDim];
    Coord_p1_Saved[iPoint] = new double[nDim];
    GridVel_Saved[iPoint] = new double[nDim];
    GridVel_Grad_Saved[iPoint] = new double*[nDim];
    for (int iDim = 0; iDim < nDim; iDim++) {
      GridVel_Grad_Saved[iPoint][iDim] = new double[nDim];
    }
    solution_Saved[iPoint] = new double[nVar];
    solution_time_n_Saved[iPoint] = new double[nVar];
    solution_time_n1_Saved[iPoint] = new double[nVar];
  }
}

Precice::~Precice(void) {
  for (int i = 0; i < localNumberWetSurfaces; i++) {
    delete [] vertexIDs[i];
  }
  delete [] vertexIDs;
  delete [] displID;
  delete [] forceID;
  delete [] displDeltaID;
  delete [] forces;
  delete [] displacementDeltas;
  delete [] valueMarkerWet;
  delete [] vertexSize;
  delete [] indexMarkerWetMappingLocalToGlobal;

  for (int iPoint = 0; iPoint < nPoint; iPoint++) {
    delete [] Coord_Saved[iPoint];
    delete [] Coord_n_Saved[iPoint];
    delete [] Coord_n1_Saved[iPoint];
    delete [] Coord_p1_Saved[iPoint];
    delete [] GridVel_Saved[iPoint];
    for (int iDim = 0; iDim < nDim; iDim++) {
      delete [] GridVel_Grad_Saved[iPoint][iDim];
    }
    delete [] GridVel_Grad_Saved[iPoint];
    delete [] solution_Saved[iPoint];
    delete [] solution_time_n_Saved[iPoint];
    delete [] solution_time_n1_Saved[iPoint];
  }
  delete [] Coord_Saved;
  delete [] Coord_n_Saved;
  delete [] Coord_n1_Saved;
  delete [] Coord_p1_Saved;
  delete [] GridVel_Saved;
  delete [] GridVel_Grad_Saved;
  delete [] solution_Saved;
  delete [] solution_time_n_Saved;
  delete [] solution_time_n1_Saved;
}


void Precice::configure( const string& preciceConfigurationFileName ){
  if (verbosityLevel_high) {
    cout << "Process #" << solverProcessIndex << "/" << solverProcessSize-1 << ": Configuring preCICE..." << endl;
  }
  solverInterface.configure( preciceConfigurationFileName );
  if (verbosityLevel_high) {
    cout << "Process #" << solverProcessIndex << "/" << solverProcessSize-1 << ": ...done configuring preCICE!" << endl;
  }
}

double Precice::initialize(){
  if (verbosityLevel_high) {
    cout << "Process #" << solverProcessIndex << "/" << solverProcessSize-1 << ": Initializing preCICE..." << endl;
  }

  //Checking for dimensional consistency of SU2 and preCICE - Exit if not consistent
  if(solverInterface.getDimensions() != geometry_container[ZONE_0][MESH_0]->GetnDim()){
    cout << "Dimensions of SU2 and preCICE are not equal! Now exiting..." << endl;
    exit(EXIT_FAILURE);
  }
  int *meshID;
  //Checking for number of wet surfaces - Exit if not cat least one wet surface defined
  if(globalNumberWetSurfaces < 1){
    cout << "There must be at least one wet surface! Now exiting..." << endl;
    exit(EXIT_FAILURE);
  } else {
    meshID = new int[globalNumberWetSurfaces];
    forceID = new int[globalNumberWetSurfaces];
    displDeltaID = new int[globalNumberWetSurfaces];
    for (int i = 0; i < globalNumberWetSurfaces; i++) {
      //Get preCICE meshIDs
      meshID[i] = solverInterface.getMeshID("SU2_Mesh" + to_string(i));
    }
  }

  //Determine the number of wet surfaces, that this process is working on, then loop over this number for all respective preCICE-related tasks
  for (int i = 0; i < globalNumberWetSurfaces; i++) {
    if (config_container[ZONE_0]->GetMarker_All_TagBound(config_container[ZONE_0]->GetpreCICE_WetSurfaceMarkerName() + to_string(i)) == -1) {
      cout << "Process #" << solverProcessIndex << "/" << solverProcessSize-1 << ": Does not work on " << config_container[ZONE_0]->GetpreCICE_WetSurfaceMarkerName() << i << endl;
    } else {
      localNumberWetSurfaces++;
    }
  }
  if (localNumberWetSurfaces < 1) {
    cout << "Process #" << solverProcessIndex << "/" << solverProcessSize-1 << ": Does not work on the wet surface at all." << endl;
    processWorkingOnWetSurface = false;
  }

  if (processWorkingOnWetSurface) {
    //Store the wet surface marker values in an array, which has the size equal to the number of wet surfaces actually being worked on by this process
    valueMarkerWet = new short[localNumberWetSurfaces];
    indexMarkerWetMappingLocalToGlobal = new short[localNumberWetSurfaces];
    int j = 0;
    for (int i = 0; i < globalNumberWetSurfaces; i++) {
      if (config_container[ZONE_0]->GetMarker_All_TagBound(config_container[ZONE_0]->GetpreCICE_WetSurfaceMarkerName() + to_string(i)) != -1) {
        valueMarkerWet[j] = config_container[ZONE_0]->GetMarker_All_TagBound(config_container[ZONE_0]->GetpreCICE_WetSurfaceMarkerName() + to_string(i));
        indexMarkerWetMappingLocalToGlobal[j] = i;
        j++;
      }
    }
    vertexIDs = new int*[localNumberWetSurfaces];
  }

  if (processWorkingOnWetSurface) {
    vertexSize = new unsigned long[localNumberWetSurfaces];
    for (int i = 0; i < localNumberWetSurfaces; i++) {
      vertexSize[i] = geometry_container[ZONE_0][MESH_0]->nVertex[valueMarkerWet[i]];

      double coupleNodeCoord[vertexSize[i]][nDim];  /*--- coordinates of all nodes at the wet surface ---*/

      unsigned long iNode;  /*--- variable for storing the node indices - one at the time ---*/
      //Loop over the vertices of the (each) boundary
      for (int iVertex = 0; iVertex < vertexSize[i]; iVertex++) {
        //Get node number (= index) to vertex (= node)
        iNode = geometry_container[ZONE_0][MESH_0]->vertex[valueMarkerWet[i]][iVertex]->GetNode();

        //Get coordinates for nodes
        for (int iDim = 0; iDim < nDim; iDim++) {
          coupleNodeCoord[iVertex][iDim] = geometry_container[ZONE_0][MESH_0]->node[iNode]->GetCoord(iDim);
          if (verbosityLevel_high) {
            cout << "Process #" << solverProcessIndex << "/" << solverProcessSize-1 << ": Initial coordinates of node (local index, global index, node color): (" << iVertex << ", " << iNode << ", " << geometry_container[ZONE_0][MESH_0]->node[iNode]->GetColor() << "): " << coupleNodeCoord[iVertex][iDim] << endl; /*--- for debugging purposes ---*/
          }
        }
      }
      //preCICE conform the coordinates of vertices (= points = nodes) at wet surface
      double coords[vertexSize[i]*nDim];
      for (int iVertex = 0; iVertex < vertexSize[i]; iVertex++) {
        for (int iDim = 0; iDim < nDim; iDim++) {
          coords[iVertex*nDim + iDim] = coupleNodeCoord[iVertex][iDim];
        }
      }

      //preCICE internal
      vertexIDs[i] = new int[vertexSize[i]];

      solverInterface.setMeshVertices(meshID[indexMarkerWetMappingLocalToGlobal[i]], vertexSize[i], coords, vertexIDs[i]);
      forceID[indexMarkerWetMappingLocalToGlobal[i]] = solverInterface.getDataID("Forces" + to_string(indexMarkerWetMappingLocalToGlobal[i]), meshID[indexMarkerWetMappingLocalToGlobal[i]]);
      displDeltaID[indexMarkerWetMappingLocalToGlobal[i]] = solverInterface.getDataID("DisplacementDeltas" + to_string(indexMarkerWetMappingLocalToGlobal[i]), meshID[indexMarkerWetMappingLocalToGlobal[i]]);
    }
    for (int i = 0; i < globalNumberWetSurfaces; i++) {
      bool flag = false;
      for (int j = 0; j < localNumberWetSurfaces; j++) {
        if (indexMarkerWetMappingLocalToGlobal[j] == i) {
          flag = true;
        }
      }
      if (!flag) {
        solverInterface.setMeshVertices(meshID[i], 0, NULL, NULL);
        forceID[i] = solverInterface.getDataID("Forces" + to_string(i), meshID[i]);
        displDeltaID[i] = solverInterface.getDataID("DisplacementDeltas" + to_string(i), meshID[i]);
      }
    }
  } else {
    for (int i = 0; i < globalNumberWetSurfaces; i++) {
      solverInterface.setMeshVertices(meshID[i], 0, NULL, NULL);
      forceID[i] = solverInterface.getDataID("Forces" + to_string(i), meshID[i]);
      displDeltaID[i] = solverInterface.getDataID("DisplacementDeltas" + to_string(i), meshID[i]);
    }
  }

  if (verbosityLevel_high) {
    cout << "Process #" << solverProcessIndex << "/" << solverProcessSize-1 << ": There is grid movement (expected: 1): " << config_container[ZONE_0]->GetGrid_Movement() << endl; /*--- for debugging purposes ---*/
    cout << "Process #" << solverProcessIndex << "/" << solverProcessSize-1 << ": Kind of grid movement (expected: 13): " << config_container[ZONE_0]->GetKind_GridMovement(ZONE_0) << endl;  /*--- for debugging purposes ---*/
  }

  double precice_dt;  /*--- preCICE timestep size ---*/
  precice_dt = solverInterface.initialize();
  if (verbosityLevel_high) {
    cout << "Process #" << solverProcessIndex << "/" << solverProcessSize-1 << ": ...done initializing preCICE!" << endl;
  }
  delete [] meshID;
  return precice_dt;
}

double Precice::advance( double computedTimestepLength ){
  if (processWorkingOnWetSurface) {
    if (verbosityLevel_high) {
      cout << "Process #" << solverProcessIndex << "/" << solverProcessSize-1 << ": Advancing preCICE..." << endl;
    }

    //Get physical simulation information
    bool incompressible = (config_container[ZONE_0]->GetKind_Regime() == INCOMPRESSIBLE);
    bool viscous_flow = ((config_container[ZONE_0]->GetKind_Solver() == NAVIER_STOKES) || (config_container[ZONE_0]->GetKind_Solver() == RANS));

    //Compute factorForces for redimensionalizing forces ("ND" = Non-Dimensional)
    double* Velocity_Real = config_container[ZONE_0]->GetVelocity_FreeStream();
    double Density_Real = config_container[ZONE_0]->GetDensity_FreeStream();
    double* Velocity_ND = config_container[ZONE_0]->GetVelocity_FreeStreamND();
    double Density_ND = config_container[ZONE_0]->GetDensity_FreeStreamND();
    double Velocity2_Real = 0.0;  /*--- denotes squared real velocity ---*/
    double Velocity2_ND = 0.0;  /*--- denotes squared non-dimensional velocity ---*/
    //Compute squared values
    for (int iDim = 0; iDim < nDim; iDim++){
      Velocity2_Real += Velocity_Real[iDim]*Velocity_Real[iDim];
      Velocity2_ND += Velocity_ND[iDim]*Velocity_ND[iDim];
    }
    //Compute factor for redimensionalizing forces
    double factorForces = Density_Real*Velocity2_Real/(Density_ND*Velocity2_ND);
    if (verbosityLevel_high) {
      cout << "Process #" << solverProcessIndex << "/" << solverProcessSize-1 << ": Factor for (non-/re-)dimensionalization of forces: " << factorForces << endl;  /*--- for debugging purposes ---*/
    }

    for (int i = 0; i < localNumberWetSurfaces; i++) {
      if (verbosityLevel_high) {
        //1. Compute forces
        cout << "Process #" << solverProcessIndex << "/" << solverProcessSize-1 << ": Advancing preCICE: Computing forces for " << config_container[ZONE_0]->GetpreCICE_WetSurfaceMarkerName() << indexMarkerWetMappingLocalToGlobal[i] << "..." << endl;
      }
      //Some variables to be used:
      unsigned long nodeVertex[vertexSize[i]];
      double normalsVertex[vertexSize[i]][nDim];
      double normalsVertex_Unit[vertexSize[i]][nDim];
      double Area;
      double Pn = 0.0;  /*--- denotes pressure at a node ---*/
      double Pinf = 0.0;  /*--- denotes environmental (farfield) pressure ---*/
      double** Grad_PrimVar = NULL; /*--- denotes (u.A. velocity) gradients needed for computation of viscous forces ---*/
      double Viscosity = 0.0;
      double Tau[3][3];
      double TauElem[3];
      double forces_su2[vertexSize[i]][nDim];  /*--- forces will be stored such, before converting to simple array ---*/

      /*--- Loop over vertices of coupled boundary ---*/
      for (int iVertex = 0; iVertex < vertexSize[i]; iVertex++) {
        //Get node number (= index) to vertex (= node)
        nodeVertex[iVertex] = geometry_container[ZONE_0][MESH_0]->vertex[valueMarkerWet[i]][iVertex]->GetNode(); /*--- Store all nodes (indices) in a vector ---*/
        // Get normal vector
        for (int iDim = 0; iDim < nDim; iDim++){
          normalsVertex[iVertex][iDim] = (geometry_container[ZONE_0][MESH_0]->vertex[valueMarkerWet[i]][iVertex]->GetNormal())[iDim];
        }
        // Unit normals
        Area = 0.0;
        for (int iDim = 0; iDim < nDim; iDim++) {
          Area += normalsVertex[iVertex][iDim]*normalsVertex[iVertex][iDim];
        }
        Area = sqrt(Area);
        for (int iDim = 0; iDim < nDim; iDim++) {
          normalsVertex_Unit[iVertex][iDim] = normalsVertex[iVertex][iDim]/Area;
        }
        // Get the values of pressure and viscosity - depending on the kind of problem
        if (incompressible){
          Pn = solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[nodeVertex[iVertex]]->GetPressureInc();
          Pinf = solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetPressure_Inf();
          if (viscous_flow){
            Grad_PrimVar = solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[nodeVertex[iVertex]]->GetGradient_Primitive();
            Viscosity = solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[nodeVertex[iVertex]]->GetLaminarViscosityInc();
          }
        }
        else {
          Pinf = solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetPressure_Inf();
          Pn = solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[nodeVertex[iVertex]]->GetPressure();
          if (viscous_flow){
            Grad_PrimVar = solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[nodeVertex[iVertex]]->GetGradient_Primitive();
            Viscosity = solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[nodeVertex[iVertex]]->GetLaminarViscosity();
          }
        }

        // Calculate the forces_su2 in the nodes for the inviscid term --> Units of force (non-dimensional).
        for (int iDim = 0; iDim < nDim; iDim++) {
          forces_su2[iVertex][iDim] = -(Pn-Pinf)*normalsVertex[iVertex][iDim];
        }
        // Calculate the forces_su2 in the nodes for the viscous term
        if (viscous_flow){
          // Divergence of the velocity
          double div_vel = 0.0;
          for (int iDim = 0; iDim < nDim; iDim++){
            div_vel += Grad_PrimVar[iDim+1][iDim];
          }
          if (incompressible){
            div_vel = 0.0;  /*--- incompressible flow is divergence-free ---*/
          }
          for (int iDim = 0; iDim < nDim; iDim++) {
            for (int jDim = 0 ; jDim < nDim; jDim++) {
              // Dirac delta
              double Delta = 0.0;
              if (iDim == jDim){
                Delta = 1.0;
              }
              // Viscous stress
              Tau[iDim][jDim] = Viscosity*(Grad_PrimVar[jDim+1][iDim] + Grad_PrimVar[iDim+1][jDim]) -
              2/3*Viscosity*div_vel*Delta;
              // Add Viscous component in the forces_su2 vector --> Units of force (non-dimensional).
              forces_su2[iVertex][iDim] += Tau[iDim][jDim]*normalsVertex[iVertex][jDim];
            }
          }
        }
        // Rescale forces_su2 to SI units
        for (int iDim = 0; iDim < nDim; iDim++) {
          forces_su2[iVertex][iDim] = forces_su2[iVertex][iDim]*factorForces;
        }
      }
      //convert forces_su2 into forces
      forces = new double[vertexSize[i]*nDim];
      for (int iVertex = 0; iVertex < vertexSize[i]; iVertex++) {
        for (int iDim = 0; iDim < nDim; iDim++) {
          //Do not write forces for duplicate nodes! -> Check wether the color of the node matches the MPI-rank of this process. Only write forces, if node originally belongs to this process.
          if (geometry_container[ZONE_0][MESH_0]->node[nodeVertex[iVertex]]->GetColor() == solverProcessIndex) {
            forces[iVertex*nDim + iDim] = forces_su2[iVertex][iDim];
          }
          else{
            forces[iVertex*nDim + iDim] = 0;
          }
        }
      }
      if (verbosityLevel_high) {
        cout << "Process #" << solverProcessIndex << "/" << solverProcessSize-1 << ": Advancing preCICE: ...done computing forces for " << config_container[ZONE_0]->GetpreCICE_WetSurfaceMarkerName() << indexMarkerWetMappingLocalToGlobal[i] << endl;
      }

      //2. Write forces
      if (verbosityLevel_high) {
        cout << "Process #" << solverProcessIndex << "/" << solverProcessSize-1 << ": Advancing preCICE: Writing forces for " << config_container[ZONE_0]->GetpreCICE_WetSurfaceMarkerName() << indexMarkerWetMappingLocalToGlobal[i] << "..." << endl;
      }
      //Load Ramping functionality: Reduce force vector before transferring by a ramping factor, which increases with the number of elapsed time steps; Achtung: ExtIter beginnt bei 0 (ohne Restart) und bei einem Restart (Startlösung) nicht bei 0, sondern bei der Startiterationsnummer
      if (config_container[ZONE_0]->GetpreCICE_LoadRamping() && ((config_container[ZONE_0]->GetExtIter() - config_container[ZONE_0]->GetUnst_RestartIter()) < config_container[ZONE_0]->GetpreCICE_LoadRampingDuration())) {
        if (verbosityLevel_high) {
          cout << "Process #" << solverProcessIndex << "/" << solverProcessSize-1 << ": Load ramping factor in preCICE: " << config_container[ZONE_0]->GetExtIter() - config_container[ZONE_0]->GetUnst_RestartIter() + 1 << "/" << config_container[ZONE_0]->GetpreCICE_LoadRampingDuration() << endl;
        }
        *forces = *forces * ((config_container[ZONE_0]->GetExtIter() - config_container[ZONE_0]->GetUnst_RestartIter()) + 1) / config_container[ZONE_0]->GetpreCICE_LoadRampingDuration();
      }
      solverInterface.writeBlockVectorData(forceID[indexMarkerWetMappingLocalToGlobal[i]], vertexSize[i], vertexIDs[i], forces);
      if (verbosityLevel_high) {
        cout << "Process #" << solverProcessIndex << "/" << solverProcessSize-1 << ": Advancing preCICE: ...done writing forces for " << config_container[ZONE_0]->GetpreCICE_WetSurfaceMarkerName() << indexMarkerWetMappingLocalToGlobal[i] << "." << endl;
      }
      if (forces != NULL){
        delete [] forces;
      }
    }

    //3. Advance solverInterface
    if (verbosityLevel_high) {
      cout << "Process #" << solverProcessIndex << "/" << solverProcessSize-1 << ": Advancing preCICE: Advancing SolverInterface..." << endl;
    }
    double max_precice_dt;
    max_precice_dt = solverInterface.advance( computedTimestepLength );
    if (verbosityLevel_high) {
      cout << "Process #" << solverProcessIndex << "/" << solverProcessSize-1 << ": Advancing preCICE: ...done advancing SolverInterface." << endl;
    }

    // displacements = new double[vertexSize*nDim]; //TODO: Delete later
    for (int i = 0; i < localNumberWetSurfaces; i++) {
      //4. Read displacements/displacementDeltas
      if (verbosityLevel_high) {
        cout << "Process #" << solverProcessIndex << "/" << solverProcessSize-1 << ": Advancing preCICE: Reading displacement deltas for " << config_container[ZONE_0]->GetpreCICE_WetSurfaceMarkerName() << indexMarkerWetMappingLocalToGlobal[i] << "..." << endl;
      }
      double displacementDeltas_su2[vertexSize[i]][nDim]; /*--- displacementDeltas will be stored such, before converting to simple array ---*/
      displacementDeltas = new double[vertexSize[i]*nDim];
      solverInterface.readBlockVectorData(displDeltaID[indexMarkerWetMappingLocalToGlobal[i]], vertexSize[i], vertexIDs[i], displacementDeltas);
      if (verbosityLevel_high) {
        cout << "Process #" << solverProcessIndex << "/" << solverProcessSize-1 << ": Advancing preCICE: ...done reading displacement deltas for " << config_container[ZONE_0]->GetpreCICE_WetSurfaceMarkerName() << indexMarkerWetMappingLocalToGlobal[i] << "." << endl;
      }

      //5. Set displacements/displacementDeltas
      if (verbosityLevel_high) {
        cout << "Process #" << solverProcessIndex << "/" << solverProcessSize-1 << ": Advancing preCICE: Setting displacement deltas for " << config_container[ZONE_0]->GetpreCICE_WetSurfaceMarkerName() << indexMarkerWetMappingLocalToGlobal[i] << "..." << endl;
      }
      //convert displacementDeltas into displacementDeltas_su2
      for (int iVertex = 0; iVertex < vertexSize[i]; iVertex++) {
        for (int iDim = 0; iDim < nDim; iDim++) {
          displacementDeltas_su2[iVertex][iDim] = displacementDeltas[iVertex*nDim + iDim];
        }
      }
      if (displacementDeltas != NULL) {
        delete [] displacementDeltas;
      }

      //Set change of coordinates (i.e. displacementDeltas)
      for (int iVertex = 0; iVertex < vertexSize[i]; iVertex++) {
        geometry_container[ZONE_0][MESH_0]->vertex[valueMarkerWet[i]][iVertex]->SetVarCoord(displacementDeltas_su2[iVertex]);
      }
      if (verbosityLevel_high) {
        cout << "Process #" << solverProcessIndex << "/" << solverProcessSize-1 << ": Advancing preCICE: ...done setting displacement deltas for " << config_container[ZONE_0]->GetpreCICE_WetSurfaceMarkerName() << indexMarkerWetMappingLocalToGlobal[i] << "." << endl;
      }
    }

    if (verbosityLevel_high) {
      cout << "Process #" << solverProcessIndex << "/" << solverProcessSize-1 << ": ...done advancing preCICE!" << endl;
    }
    return max_precice_dt;
  }
  else{
    if (verbosityLevel_high) {
      cout << "Process #" << solverProcessIndex << "/" << solverProcessSize-1 << ": Advancing preCICE..." << endl;
    }
    //3. Advance solverInterface
    if (verbosityLevel_high) {
      cout << "Process #" << solverProcessIndex << "/" << solverProcessSize-1 << ": Advancing preCICE: Advancing SolverInterface..." << endl;
    }
    double max_precice_dt;
    max_precice_dt = solverInterface.advance( computedTimestepLength );
    if (verbosityLevel_high) {
      cout << "Process #" << solverProcessIndex << "/" << solverProcessSize-1 << ": Advancing preCICE: ...done advancing SolverInterface." << endl;
      cout << "Process #" << solverProcessIndex << "/" << solverProcessSize-1 << ": ...done advancing preCICE!" << endl;
    }
    return max_precice_dt;
  }
}

bool Precice::isCouplingOngoing(){
  return solverInterface.isCouplingOngoing();
}

bool Precice::isActionRequired( const string& action ){
  return solverInterface.isActionRequired(action);
}

const string& Precice::getCowic(){
  return cowic;
}

const string& Precice::getCoric(){
  return coric;
}

void Precice::saveOldState( bool *StopCalc, double *dt ){

  for (int iPoint = 0; iPoint < nPoint; iPoint++) {
    for (int iVar = 0; iVar < nVar; iVar++) {
      //Save solutions at last and current time step
      solution_Saved[iPoint][iVar] = (solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->GetSolution())[iVar];
      solution_time_n_Saved[iPoint][iVar] = (solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->GetSolution_time_n())[iVar];
      solution_time_n1_Saved[iPoint][iVar] = (solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->GetSolution_time_n1())[iVar];
    }
    for (int iDim = 0; iDim < nDim; iDim++) {
      //Save coordinates at last, current and next time step
      Coord_Saved[iPoint][iDim] = (geometry_container[ZONE_0][MESH_0]->node[iPoint]->GetCoord())[iDim];
      Coord_n_Saved[iPoint][iDim] = (geometry_container[ZONE_0][MESH_0]->node[iPoint]->GetCoord_n())[iDim];
      Coord_n1_Saved[iPoint][iDim] = (geometry_container[ZONE_0][MESH_0]->node[iPoint]->GetCoord_n1())[iDim];
      Coord_p1_Saved[iPoint][iDim] = (geometry_container[ZONE_0][MESH_0]->node[iPoint]->GetCoord_p1())[iDim];
      //Save grid velocity
      GridVel_Saved[iPoint][iDim] = (geometry_container[ZONE_0][MESH_0]->node[iPoint]->GetGridVel())[iDim];
      for (int jDim = 0; jDim < nDim; jDim++) {
        //Save grid velocity gradient
        GridVel_Grad_Saved[iPoint][iDim][jDim] = (geometry_container[ZONE_0][MESH_0]->node[iPoint]->GetGridVel_Grad())[iDim][jDim];
      }
    }
  }

  //Save wether simulation should be stopped after the current iteration
  StopCalc_savedState = *StopCalc;
  //Save the time step size
  dt_savedState = *dt;
  //Writing task has been fulfilled successfully
  solverInterface.fulfilledAction(cowic);
}

void Precice::reloadOldState(bool *StopCalc, double *dt){
  for (int iPoint = 0; iPoint < nPoint; iPoint++) {
    //Reload solutions at last and current time step
    solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->SetSolution(solution_Saved[iPoint]);
    solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->Set_Solution_time_n(solution_time_n_Saved[iPoint]);
    solver_container[ZONE_0][MESH_0][FLOW_SOL]->node[iPoint]->Set_Solution_time_n1(solution_time_n1_Saved[iPoint]);

    //Reload coordinates at last, current and next time step
    geometry_container[ZONE_0][MESH_0]->node[iPoint]->SetCoord(Coord_n1_Saved[iPoint]);
    geometry_container[ZONE_0][MESH_0]->node[iPoint]->SetCoord_n();
    geometry_container[ZONE_0][MESH_0]->node[iPoint]->SetCoord_n1();
    geometry_container[ZONE_0][MESH_0]->node[iPoint]->SetCoord(Coord_n_Saved[iPoint]);
    geometry_container[ZONE_0][MESH_0]->node[iPoint]->SetCoord_n();
    geometry_container[ZONE_0][MESH_0]->node[iPoint]->SetCoord_p1(Coord_p1_Saved[iPoint]);
    geometry_container[ZONE_0][MESH_0]->node[iPoint]->SetCoord(Coord_Saved[iPoint]);

    //Reload grid velocity
    geometry_container[ZONE_0][MESH_0]->node[iPoint]->SetGridVel(GridVel_Saved[iPoint]);

    //Reload grid velocity gradient
    for (int iDim = 0; iDim < nDim; iDim++) {
      for (int jDim = 0; jDim < nDim; jDim++) {
        geometry_container[ZONE_0][MESH_0]->node[iPoint]->SetGridVel_Grad(iDim, jDim, GridVel_Grad_Saved[iPoint][iDim][jDim]);
      }
    }
  }

  //Reload wether simulation should be stopped after current iteration
  *StopCalc = StopCalc_savedState;
  //Reload the time step size
  *dt = dt_savedState;
  //Reading task has been fulfilled successfully
  solverInterface.fulfilledAction(coric);
}

void Precice::finalize(){
  cout << "Process #" << solverProcessIndex << "/" << solverProcessSize-1 << ": Finalizing preCICE..." << endl;
  solverInterface.finalize();
  cout << "Process #" << solverProcessIndex << "/" << solverProcessSize-1 << ": Done finalizing preCICE!" << endl;
}
