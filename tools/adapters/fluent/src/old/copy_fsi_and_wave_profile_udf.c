#include "udf.h"

#ifdef __cplusplus
extern "C" {
#endif
#include "wave_profile_udf.h"
#ifdef __cplusplus
}
#endif

#include "dynamesh_tools.h"
#include "../../../../src/precice/adapters/c/SolverInterfaceC.h"
#include "../../../../src/precice/adapters/c/Constants.h"
#include <float.h>
#include <math.h>
#include <string.h>
#include <stdio.h>

#ifdef BOOL_TRUE
#error BOOL_TRUE already defined!
#endif

#ifdef BOOL_FALSE
#error BOOL_TRUE already defined!
#endif

#define BOOL_TRUE  1
#define BOOL_FALSE 0

#if ND_ND == 2
#define norm(a, b) sqrt((a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1]))
#else
#error Not implemented!
#endif

/* Wave variables */
const real PI = 3.1415926535897932;
const real GRAVITY = 9.81;
const real SHIFT = 3.1415926535897932*0.5;
int counter = 0;
real wave_amplitude = 0.0;
real omega = 0.0;
real k = 0.0;
real surface_level = 0.0;

/* FSI variables */
double timestep_limit = 0.0;
char* const createIterationCheckpoint = "write-iteration-checkpoint";
char* const readIterationCheckpoint = "read-iteration-checkpoint";
double* forces = NULL;
int* force_indices = NULL;
int skip_grid_motion = BOOL_TRUE;
int gather_write_positions = BOOL_TRUE;
int gather_read_positions = BOOL_TRUE;
int thread_index = 0;
int dynamic_thread_size = 0;
int wet_edge_size = 0;
int wet_nodes_size = 0;
int boundary_nodes_size = 0;
int deformable_nodes_size = 0;
int moved_nodes_counter = 0;
double* initial_coords = NULL;
double* boundary_coords = NULL; /* MESH MOVEMENT */
double* displacements = NULL;
int* displ_indices = NULL;
int* dynamic_thread_node_size = NULL;
double* c_matrix = NULL;
double* x_coeff_vector = NULL;
double* y_coeff_vector = NULL;
double* b_vector = NULL;
int* pivots_vector = NULL;
int process_id = -1;
int comm_size = -1;
int require_create_checkpoint = BOOL_FALSE;

/**
 * Sets up initial distribution of water and air.
 */
void init_phases(Domain* domain);

/**
 * Inits FSI related stuff.
 */
void init_fsi();

/**
 * Returns the phase ratio (by water) of a boundary cell.
 */
real get_ratio_edge_is_covered (
  face_t  face,
  Thread* face_thread,
  real    wave_height );

DEFINE_INIT(init,domain)
{
  Message("\nEntering INIT\n");
  init_phases(domain);
  init_fsi();
  Message("Leaving INIT\n");
}

DEFINE_PROFILE(MP_wave_x_velocity,face_thread,var_index)
{
  /*Message ( "\nEntering DEFINE_PROFILE(wave_x_velocity)\n");*/
  real time = CURRENT_TIME;
  real centroid_coords[ND_ND];
  face_t face;
  real height = 0.0;
  real vel = 0.0;
  real wave_height = wave_amplitude * cos(-1.0 * omega * time - SHIFT);
  real vel_amplitude = (wave_amplitude * GRAVITY * k * cos(-1.0 * omega * time - SHIFT)) / (cosh(k * surface_level) * omega);

# if !RP_HOST /* Serial or node */
  begin_f_loop(face, face_thread){
    if (PRINCIPAL_FACE_P(face, face_thread)){
      F_CENTROID(centroid_coords, face, face_thread);
      height = centroid_coords[1] - surface_level;
      if (height <= wave_height){
        vel = vel_amplitude * cosh(k * (height + surface_level));
      }
      else {
        vel = vel_amplitude * cosh(k * (wave_height + surface_level));
      }
      /*Message("Setting x vel = %f at face centroid_y = %f\n", vel, height);*/
      F_PROFILE(face, face_thread, var_index) = vel;
    }
  } end_f_loop(face, face_thread)
# endif
  /*Message ( "Leaving DEFINE_PROFILE(wave_x_velocity)\n" );*/
}

DEFINE_PROFILE(MP_wave_y_velocity, face_thread, var_index)
{
  /*Message ( "\nEntering DEFINE_PROFILE(wave_y_velocity)\n" );*/
  real time = CURRENT_TIME;
  real centroid_coords[ND_ND];
  face_t face;
  real height = 0.0;
  real vel = 0.0;
  real wave_height = wave_amplitude * cos(-1.0 * omega * time - SHIFT);
  real vel_amplitude = (wave_amplitude * GRAVITY * k * sin(-1.0 * omega * time - SHIFT)) / (cosh(k * surface_level) * omega);

# if !RP_HOST /* Serial or node */
  begin_f_loop(face, face_thread){
    if (PRINCIPAL_FACE_P(face, face_thread)){
      F_CENTROID(centroid_coords, face, face_thread);
      height = centroid_coords[1] - surface_level;
      if (height <= wave_height){
        vel = vel_amplitude * sinh(k * (height + surface_level));
        /*vel = 0.2 * sin(-1.0 * omega * time  - SHIFT);*/
      }
      else {
        vel = vel_amplitude * sinh(k * (wave_height + surface_level));
      }
      /*Message("Setting y vel = %f at face centroid_y = %f\n", vel, height);*/
      F_PROFILE(face, face_thread, var_index) = vel;
    }
  } end_f_loop(face, face_thread)
# endif
  /*Message ( "Leaving DEFINE_PROFILE(wave_y_velocity)\n" );*/
}

DEFINE_PROFILE(MP_wave_volume_fraction, face_thread, var_index)
{
  /*Message ( "\nEntering DEFINE_PROFILE(wave_volume_fraction)\n" );*/
  real time = CURRENT_TIME;
  real centroid_coords[ND_ND];
  face_t face;
  real water_fraction = 0.0;
  real wave_height = 0.0;
  real centroid_y = 0.0;
  real upper_y = 0.0;
  real lower_y = 0.0;
  real area_vector[ND_ND];
  real area = 0.0;
  wave_height = wave_amplitude * cos(-1.0 * omega * time - SHIFT);

# if !RP_HOST /* Serial or node */
  begin_f_loop(face, face_thread){
    if (PRINCIPAL_FACE_P(face, face_thread)){
      F_CENTROID(centroid_coords, face, face_thread);
      F_AREA(area_vector, face, face_thread);
      area = NV_MAG(area_vector);
      centroid_y = centroid_coords[1];
      upper_y = centroid_y + area/2.0;
      lower_y = centroid_y - area/2.0;
      if (lower_y > surface_level + wave_height){
         water_fraction = 0.0;
      }
      else if (upper_y < surface_level + wave_height){
         water_fraction = 1.0;
      }
      else {
         water_fraction = (surface_level + wave_height - lower_y) / area;
      }
      /*water_fraction = height < surface_level + wave_height ? 1.0 : 0.0;*/
      /*Message("Setting water fraction = %f at face centroid_y = %f\n", water_fraction, centroid_y);*/
      F_PROFILE(face, face_thread, var_index) = water_fraction;
    }
  } end_f_loop(face, face_thread)
# endif
  /*Message ( "Leaving DEFINE_PROFILE(wave_volume_fraction)\n" );*/
}

DEFINE_ON_DEMAND(write_and_advance)
{
  Message("(%d) Entering ON_DEMAND(write_and_andvance)\n", process_id);
  int meshID = precice_getMeshID("WetSurface");
  int forceID = precice_getDataID("Forces");
  real pressure_force[ND_ND];
  real viscous_force[ND_ND];
  double total_force[ND_ND];
  double max_force = 0.0;
  int i, j; /* Loop counter */
  double center[ND_ND];
  real area[ND_ND];
  int ongoing;
  int subcycling = ! precice_isWriteDataRequired(CURRENT_TIMESTEP);
  Dynamic_Thread* dynamic_thread = NULL;
  int thread_counter = 0;

# if !RP_HOST /* Serial or Node */
  if (! subcycling){
    /*Message ( "   Write force data\n" );*/
    Domain *domain = Get_Domain(1);
    if (domain == NULL){
      Message("      (%d) ERROR: domain == NULL\n", process_id);
      exit(1);
    }
    /*Message ( "      Domain name = %s", domain->name );*/
    if (domain->dynamic_threads == NULL){
      Message("      (%d) ERROR: domain.dynamic_threads == NULL\n", process_id);
      exit(1);
    }

    if (gather_write_positions){
      Message("      (%d) Counting wet edges...", process_id);
      dynamic_thread = domain->dynamic_threads;
      thread_counter = 0;
      while (dynamic_thread != NULL){
        if (strncmp("gridmotions", dynamic_thread->profile_udf_name, 11) == 0){
          Message("\n      (%d) Thread index %d\n", process_id, thread_counter);
          Thread* face_thread  = DT_THREAD(dynamic_thread);
          if (face_thread == NULL){
            Message("      (%d) ERROR: face_thread == NULL\n", process_id);
            exit(1);
          }
          face_t face;
          begin_f_loop (face, face_thread){
            if (PRINCIPAL_FACE_P(face,face_thread)){
              wet_edge_size++;
            }
          } end_f_loop(face, face_thread);
          thread_counter++;
        }
        dynamic_thread = dynamic_thread->next;
      }
      Message("      (%d) ...done (counted %d wet edges)", process_id, wet_edge_size);
      Message("      (%d) Allocating %d force vector values\n", process_id, wet_edge_size * ND_ND);
      if (forces != NULL){
        Message("      (%d) ERROR: Forces vector allocated multiple times!\n", process_id);
      }
      forces = (double*) malloc(wet_edge_size * ND_ND * sizeof(double));
      Message("      (%d) Allocating %d force indices\n", process_id, wet_edge_size);
      force_indices = (int*) malloc(wet_edge_size * sizeof(int));

      Message("      (%d) Setting write positions...", process_id);
      dynamic_thread = domain->dynamic_threads;
      thread_counter = 0;
      i = 0;
      while (dynamic_thread != NULL){
        if (strncmp("gridmotions", dynamic_thread->profile_udf_name, 11) == 0){
          Message("\n      (%d) Thread index %d\n", process_id, thread_counter);
          Thread* face_thread  = DT_THREAD(dynamic_thread);
          if (face_thread == NULL){
            Message("      (%d) ERROR: face_thread == NULL\n", process_id);
            exit(1);
          }
          face_t face;
          begin_f_loop (face, face_thread){
            if (PRINCIPAL_FACE_P(face,face_thread)){
              F_CENTROID(center, face, face_thread);
              force_indices[i] = precice_setWritePosition(meshID, center);
              i++;
            }
          } end_f_loop(face, face_thread);
          thread_counter++;
        }
        dynamic_thread = dynamic_thread->next;
      }
      Message("      (%d) ...done\n", process_id);
      gather_write_positions = BOOL_FALSE;
    }

    dynamic_thread = domain->dynamic_threads;
    thread_counter = 0;
    Message("      (%d) Gather forces...\n", process_id);
    i = 0;
    while (dynamic_thread != NULL){
      if (strncmp("gridmotions", dynamic_thread->profile_udf_name, 11) == 0){
        Message("\n      (%d) Thread index %d\n", process_id, thread_counter);
        Thread* face_thread  = DT_THREAD(dynamic_thread);
        if (face_thread == NULL){
          Message("      (%d) ERROR: face_thread == NULL\n", process_id);
          exit(1);
        }
        face_t face;
        begin_f_loop (face, face_thread){
          if (PRINCIPAL_FACE_P(face,face_thread)){
            /*F_CENTROID(center, face, face_thread);*/
            F_AREA(area, face, face_thread);
            NV_VS(viscous_force, =, F_STORAGE_R_N3V(face,face_thread,SV_WALL_SHEAR),*,-1.0);
            NV_VS(pressure_force, =, area, *, F_P(face,face_thread));
            NV_VV(total_force, =, viscous_force, +, pressure_force);
            Message ( "      (%d) Area = %f, %f; Visc.Force = %f, %f; PressureForce = %f, %f; Total = %f, %f\n",
                      process_id, area[0], area[1], viscous_force[0], viscous_force[1],
                      pressure_force[0], pressure_force[1], total_force[0], total_force[1] );

            /*int index = precice_setWritePosition(meshID, center);*/
            /*Message ( "      After writing" );
            */
            /*precice_writeVectorData(forceID, index, total_force);*/
            for (j=0; j < ND_ND; j++){
              forces[i + j] = total_force[j];
              if (fabs(total_force[j]) > fabs(max_force)){
                max_force = total_force[j];
              }
            }
            i += ND_ND;
          }
        } end_f_loop(face, face_thread);
        thread_counter++;
      }
      dynamic_thread = dynamic_thread->next;
    }
    Message("      (%d) ...done (with %d force values)\n", process_id, i);
    Message("      (%d) Writing forces...\n", process_id);
    precice_writeBlockVectorData(forceID, wet_edge_size, force_indices, forces);
    Message("      (%d) ...done\n", process_id );
    Message("      (%d) Max force: %f\n", max_force);
    if (thread_counter != dynamic_thread_size){
      Message ( "      (%d) ERROR: Number of dynamic threads has changed to %d!\n", process_id, thread_counter );
      exit(1);
    }
  }
# endif /* !RP_HOST */

  timestep_limit = precice_advance(CURRENT_TIMESTEP);

  /* Read coupling state */
  ongoing = precice_isCouplingOngoing();
# if !RP_NODE
  RP_Set_Integer("udf/ongoing", ongoing);
# endif /* !RP_NODE */

  if (precice_isActionRequired(createIterationCheckpoint)){
#   if !RP_NODE
    RP_Set_Integer("udf/convergence", BOOL_TRUE);
#   endif /* !RP_NODE */
    precice_fulfilledAction(createIterationCheckpoint);
  }

  if (precice_isActionRequired(readIterationCheckpoint)){
#   if !RP_NODE
    RP_Set_Integer("udf/convergence", BOOL_FALSE);
#   endif /* !RP_NODE */
    precice_fulfilledAction(readIterationCheckpoint);
  }

# if !RP_NODE
  if (! precice_isCouplingOngoing()){
    RP_Set_Integer("udf/convergence", BOOL_TRUE);
  }
# endif /* !RP_NODE */

  if (precice_isActionRequired(precice_actionWriteSimulationCheckpoint())){
    precice_fulfilledAction(precice_actionWriteSimulationCheckpoint());
    require_create_checkpoint = BOOL_TRUE;
#   if !RP_NODE
    Message("   (%d) Writing simulation checkpoint required\n", process_id);
    RP_Set_Integer("udf/checkpoint", BOOL_TRUE);
#   endif /* !RP_NODE */
  }
# if !RP_NODE
  else {
    RP_Set_Integer("udf/checkpoint", BOOL_FALSE);
  }
# endif /* !RP_NODE */


  Message("(%d) Leaving ON_DEMAND(write_and_advance)\n", process_id);
}

DEFINE_GRID_MOTION(gridmotions,domain,dt,time,dtime)
{
  Message("\n(%d) Entering GRID_MOTION\n", process_id);
  if (skip_grid_motion){
    if (thread_index >= dynamic_thread_size-1){
      skip_grid_motion = BOOL_FALSE;
    }
    thread_index++;
    Message("   Skipping first round grid motion\n");
    return;
  }
  int meshID = precice_getMeshID("WetSurface");
# if !RP_HOST /* Serial or node */
  if (thread_index == dynamic_thread_size){
    Message ("   Reset thread index\n");
    thread_index = 0;
  }
  Message("   Thread index = %d\n", thread_index);
  Thread* face_thread  = DT_THREAD(dt);
  Node* node;
  face_t face;
  int n;
  int displID = precice_getDataID("Displacements");
  /*double displacement[ND_ND];*/
  double coords[ND_ND];
  double max_displ_delta = 0.0;
  /*int sizeNodes = 0;*/
  int node_index = 0;
  int array_index = 0;
  int dim = 0;
  int offset = 0;
  int i = 0, j = 0;
  int entries_per_line = 0; /* For RBF deformation field computation */
  double* x_i;
  double* x_j;
  double r=0.0; /* Radial distance to radial-basis-function */
  /*int subcycling = ! precice_isReadDataAvailable ();*/

  if (strncmp("gridmotions", dt->profile_udf_name, 11) != 0){
    Message("   ERROR: called gridmotions for invalid dynamic thread: %s\n",
            dt->profile_udf_name);
    exit(1);
  }

  if (face_thread == NULL){
    Message("   ERROR: face_thread == NULL\n");
    exit(1);
  }

  /*
  if ( subcycling ){
    Message ( "   Skipping grid motion due to subcycling\n" );
    return;
  }
  */

  if (gather_read_positions){
    /* Check that no node is marked updated yet */
    /*
    begin_f_loop ( face, face_thread ){
      f_node_loop ( face, face_thread, n ){
        node = F_NODE ( face, face_thread, n );
        if ( NODE_POS_UPDATED_P(node) ){
          Message ( "ERROR: Found already updated node before updating!" );
          exit (1);
        }
      }
    } end_f_loop ( face, face_thread );
    */

    /* Count not yet as updated (from other threads) marked nodes */
    dynamic_thread_node_size[thread_index] = 0;
    begin_f_loop(face, face_thread){
      if ( PRINCIPAL_FACE_P(face,face_thread) ){
        f_node_loop ( face, face_thread, n ){
          node = F_NODE ( face, face_thread, n );
          if ( NODE_POS_NEED_UPDATE(node) ){
            /*NODE_POS_UPDATED ( node );*/
            NODE_MARK(node) = 12345;
            /*sizeNodes++;*/
            wet_nodes_size++;
            dynamic_thread_node_size[thread_index]++;
          }
        }
      }
    } end_f_loop(face, face_thread);
    /*Message ( "Counted sizeNodes = %i\n", sizeNodes );*/

    /* Get initial coordinates and reset update marking */
    Message("      (%d) Reallocating %d initial positions ...\n", process_id, wet_nodes_size);
    initial_coords = (double*) realloc(initial_coords, wet_nodes_size * ND_ND * sizeof(double));
    displacements = (double*) realloc(displacements, wet_nodes_size * ND_ND * sizeof(double));
    displ_indices = (int*) realloc(displ_indices, wet_nodes_size * sizeof(int));
    /*index = 0;*/
    array_index = wet_nodes_size - dynamic_thread_node_size[thread_index];
    begin_f_loop (face, face_thread){
      if (PRINCIPAL_FACE_P(face,face_thread)){
        f_node_loop(face, face_thread, n){
          node = F_NODE(face, face_thread, n);
          if (NODE_MARK(node) == 12345){
            NODE_MARK(node) = 1;  /*Set node to need update*/
            for (dim=0; dim < ND_ND; dim++){
              /*Message("-6: dim=%d, array_index=%d\n", dim, array_index);*/
              coords[dim] = NODE_COORD(node)[dim];
              initial_coords[array_index*ND_ND+dim] = coords[dim];
            }
            node_index = precice_setReadPosition(meshID, coords);
            /*Message ( "initial_coords[%d] = %f, %f\n",
              array_index, initial_coords[array_index*ND_ND], initial_coords[array_index*ND_ND+1] );
            */
            displ_indices[array_index] = node_index;
            array_index++;
          }
        }
      }
    } end_f_loop(face, face_thread);
    Message("      (%d) Set %d (of %d) displacement read positions ...\n", process_id,
            array_index - wet_nodes_size + dynamic_thread_node_size[thread_index],
            dynamic_thread_node_size[thread_index]);

    if (thread_index == dynamic_thread_size - 1){
      gather_read_positions = BOOL_FALSE;
      /*compute_rbf_matrix();*/
    }
  }

  /*indexNode = 0;*/
  SET_DEFORMING_THREAD_FLAG(THREAD_T0(face_thread));
# endif /* !RP_HOST */

  precice_mapReadData(meshID); /* Collective call necessary */

# if !RP_HOST /* Serial or node */
  if (dynamic_thread_node_size[thread_index] > 0){
    Message("      (%d) Reading displacements...\n", process_id);
    offset = 0;
    for (i = 0; i < thread_index; i++){
      offset += dynamic_thread_node_size[i];
    }
    /*Message("      offset=%d, size=%d\n", offset, dynamic_thread_node_size[thread_index]);*/
    precice_readBlockVectorData(displID, dynamic_thread_node_size[thread_index],
        displ_indices + offset, displacements + ND_ND * offset);

    Message("      (%d) Setting displacements...\n", process_id);
    i = offset * ND_ND;
    begin_f_loop (face, face_thread){
      if (PRINCIPAL_FACE_P(face,face_thread)){
        f_node_loop (face, face_thread, n){
          node = F_NODE(face, face_thread, n);
          if (NODE_POS_NEED_UPDATE(node)){
            NODE_POS_UPDATED(node);
            /*NV_V ( coords, =, NODE_COORD_N(node) );*/
            /* precice_readVectorData ( displacementID, coords, displacement ); */
            /*int index = precice_setReadPosition ( meshID, coords );*/
            /*Message ( "Set read position at %f, %f\n",
                      *(initial_coords[thread_index] + indexNode),
                      *(initial_coords[thread_index] + indexNode + 1));
            */
            /*int index = precice_setReadPosition ( meshID, initial_coords[thread_index] + indexNode );
            precice_readVectorData ( displID, index, displacement );
            */
            /*Message ( "Scaled Displ. Deltas x=%f, y=%f", displacement[0], displacement[1] );*/
            Message("old coord = %f, %f", NODE_COORD_N(node)[0], NODE_COORD_N(node)[1] );
            Message("; displ = %f, %f", displacements[i], displacements[i+1]);

            /*NV_V ( NODE_COORD(node), =, NODE_COORD_N(node) );*/
            /*NV_V ( NODE_COORD(node), +=, displacement );*/

            for (dim=0; dim < ND_ND; dim++){
              NODE_COORD(node)[dim] = initial_coords[i+dim] + displacements[i+dim];
              displacements[i+dim] = NODE_COORD(node)[dim] - NODE_COORD_N(node)[dim]; /* store delta for rbf mesh motion */
              /*NODE_COORD(node)[dim] = NODE_COORD_N(node)[dim] + displacements[i + dim];*/
              if (fabs(displacements[i+dim]) > fabs(max_displ_delta)){
                max_displ_delta = displacements[i + dim];
              }
            }
            Message("; new coord = %f, %f\n", NODE_COORD(node)[0], NODE_COORD(node)[1]);
            i += ND_ND;
          }
        }
      }
    } end_f_loop (face, face_thread);
    Message("      (%d) ...done\n", process_id);

    /*
    if (thread_index == dynamic_thread_size-1){
      Message("      (%d) Moved nodes counter: %d, of %d\n", process_id,
              moved_nodes_counter, );
      compute_rbf_coefficients();
      moved_nodes_counter = 0;
    }
    */
  }
  Message("      (%d) Max displacement delta: %f\n", process_id, max_displ_delta);

  thread_index++;
# endif /* !RP_HOST */

# if !RP_NODE
  Message("      (%d) convergence=%d, iterate=%d, couplingOngoing=%d\n",
          process_id, RP_Get_Integer("udf/convergence"), RP_Get_Integer("udf/iterate"),
          precice_isCouplingOngoing());
  if (RP_Get_Integer("udf/convergence") && RP_Get_Integer("udf/iterate") && precice_isCouplingOngoing()){
    RP_Set_Integer("udf/convergence", BOOL_FALSE);
  }
# endif /* !RP_NODE */
  if (! precice_isCouplingOngoing()){
    precice_finalize();
  }

  Message("(%d) Leaving GRID_MOTION\n", process_id);
}

void init_phases(Domain* domain)
{
  Domain* subdomain;
  Thread* cell_thread; /* linked list of cells */
  cell_t cell;    /* one cell */
  int phase_domain_index;
  real cell_centroid[ND_ND]; /* n dimensional cell center coordinates */
  int flipzones = 0;
  int phase_sign = 0;
  real upper_phase = 0.0;
  real lower_phase = 0.0;

  int i;
  real coords[6];
  int node_index;
  Node* node;
  real phase;

# if !RP_NODE /* Serial or Host */
  surface_level = RP_Get_Float("udf/surfacelevel");
  flipzones = RP_Get_Integer("udf/flipzones");
  wave_amplitude = 0.5 * RP_Get_Float("udf/waveheight");
  omega = (2.0 * PI) / RP_Get_Float("udf/period");
  k = (2.0 * PI) / RP_Get_Float("udf/wavelength");
# endif
  host_to_node_real_4(surface_level, wave_amplitude, omega, k);
  host_to_node_int_1(flipzones);

  Message("surface level = %f\n", surface_level);
  Message("flipzones = %d\n", flipzones);
  Message("wave_amplitude = %f\n", wave_amplitude);
  Message("omega = %f\n", omega);
  Message("k = %f\n", k);
  
# if !RP_HOST /* serial or node */
  sub_domain_loop (subdomain, domain, phase_domain_index){ 
    if (phase_domain_index == 0){ /* primary phase */
      if (! flipzones){
        phase_sign = 1;
        upper_phase = 1.0;
        lower_phase = 0.0;
      }
      else {
        phase_sign = -1;
        upper_phase = 0.0;
        lower_phase = 1.0;
      }
    }
    else if (phase_domain_index == 1){
      if (! flipzones){
        phase_sign = -1;
        upper_phase = 0.0;
        lower_phase = 1.0;
      }
      else {
        phase_sign = 1;
        upper_phase = 1.0;
        lower_phase = 0.0;
      }
    }
    Message("Setting phase for upper=%f, lower=%f\n", upper_phase, lower_phase);
    thread_loop_c (cell_thread, subdomain){
      begin_c_loop_all (cell, cell_thread){
#       if ND_ND == 3
        C_CENTROID(cell_centroid, cell, cell_thread);
        if (cell_centroid[1] < surface_level){
          C_VOF(cell, cell_thread) = lower_phase;
        }
        else {
          C_VOF(cell, cell_thread) = upper_phase;
        }
#       endif /* ND_ND == 3 */
#       if ND_ND == 2
        i = 0;
        c_node_loop(cell, cell_thread, node_index){
          node = C_NODE(cell, cell_thread, node_index);
          coords[i++] = NODE_X(node);
          coords[i++] = NODE_Y(node);
        }
        /*if (i != 6){
          Message("Error: found cell with size nodes != 3\n");
          exit(1);
        }*/
        phase = phase_sign == 1
            ? 1.0 - get_ratio_triangle_is_covered(coords, surface_level)
            : get_ratio_triangle_is_covered(coords, surface_level);

        C_CENTROID(cell_centroid, cell, cell_thread);
        if (phase < 0.0 || phase > 1.0){
          Message("Error: Cell at %f,%f has phase = %f\n", cell_centroid[0], cell_centroid[1], phase);
          exit(-1);
        }

        /*if (cell_centroid[1] > 0.5 && cell_centroid[1] < 0.55){
          Message("At %f,%f: %f (coords = %f, %f, %f, %f, %f, %f)\n", cell_centroid[0], cell_centroid[1], phase,
                  coords[0], coords[1], coords[2], coords[3], coords[4], coords[5]);
        }*/
        C_VOF(cell, cell_thread) = phase;
#       endif /* ND_ND == 2 */
      } end_c_loop_all(cell_type, cell_thread)
    }
  }
# endif /* serial or node */
}

void init_fsi()
{
  Dynamic_Thread* dynamic_thread = NULL;
  Thread* face_thread = NULL;
  face_t face;
  int node_index;
  Node* node = NULL;
  int i = 0; /* Loop counter */
  int dim = 0; /* Dimension loop counter */

# if RP_NODE
  process_id = myid + 1;
  comm_size = compute_node_count + 1;
  Message("   Compute node: id=%d, size=%d\n", process_id, comm_size);
# endif
# if RP_HOST
  process_id = 0;
  comm_size = compute_node_count + 1;
  Message("   Host: id=%d, size=%d\n", process_id, comm_size);
# endif
# if !PARALLEL
  process_id = 0;
  comm_size = 1;
  Message("   (%d) Serial node: id=%d, size=%d\n", process_id, process_id, comm_size);
# endif

  Message("   (%d) Creating coupling interface\n", process_id);
  precice_createSolverInterface("Fluent", precice_nameConfiguration(),
                                process_id, comm_size);

  Message("   (%d) Initializing coupled simulation\n", process_id);
  timestep_limit = precice_initialize();
  Message("   (%d) ... done\n", process_id);

  if (precice_isActionRequired(precice_actionReadSimulationCheckpoint())){
#   if !RP_NODE
    Message("   (%d) Reading simulation checkpoint required\n", process_id);
    RP_Set_Integer("udf/checkpoint", BOOL_TRUE);
    gather_write_positions = BOOL_FALSE;
    gather_read_positions = BOOL_FALSE;
#   endif /* !RP_NODE */
    /*precice_fulfilledAction(precice_actionReadSimulationCheckpoint());*/
  }

  if (precice_isActionRequired(createIterationCheckpoint)){
    Message("   (%d) Implicit coupling\n", process_id);
#   if !RP_NODE
    RP_Set_Integer("udf/convergence", BOOL_FALSE);
    RP_Set_Integer("udf/iterate", BOOL_TRUE);
#   endif /* !RP_NODE */
    precice_fulfilledAction(createIterationCheckpoint);
  }
  else {
    Message("   (%d) Explicit coupling\n", process_id);
  }

# if !RP_HOST
  /* Count dynamic threads */
  Message ( "   (%d) counting dynamic threads: ", process_id );
  Domain *domain = Get_Domain(1);
  if (domain == NULL){
    Message("      (%d) ERROR: domain == NULL\n", process_id);
    exit(1);
  }
  if (domain->dynamic_threads == NULL){
    Message ( "      (%d) ERROR: domain.dynamic_threads == NULL\n", process_id );
    exit (1);
  }
  dynamic_thread = domain->dynamic_threads;
  while (dynamic_thread != NULL){
    if (strncmp("gridmotions", dynamic_thread->profile_udf_name, 11) == 0){
      face_thread  = DT_THREAD(dynamic_thread);
      begin_f_loop (face, face_thread){ /* Thread face loop */
        f_node_loop (face, face_thread, node_index){ /* Face node loop */
          node = F_NODE(face, face_thread, node_index);
          /*Message("   (%d) dynamic thread node mark = %d\n", process_id, NODE_MARK(node));*/
          NODE_MARK(node) = 11111;
        }
      } end_f_loop(face, face_thread)
      dynamic_thread_size++;
    }
    dynamic_thread = dynamic_thread->next;
  }
  /*displacements = (double**) malloc(dynamic_thread_size * sizeof(double*));*/
  /*displ_indices = (int**) malloc(dynamic_thread_size * sizeof(int*));*/
  dynamic_thread_node_size = (int*) malloc(dynamic_thread_size * sizeof(int));
  Message("%d\n", dynamic_thread_size);

  dynamic_thread = domain->dynamic_threads;
  while (dynamic_thread != NULL){
    if (strncmp("gridmotions", dynamic_thread->profile_udf_name, 11) == 0){
      face_thread  = DT_THREAD(dynamic_thread);
      begin_f_loop (face, face_thread){ /* Thread face loop */
        f_node_loop (face, face_thread, node_index){ /* Face node loop */
          node = F_NODE(face, face_thread, node_index);
          /*Message("   (%d) reset dynamic thread node mark = %d\n", process_id, NODE_MARK(node));*/
          NODE_MARK(node) = 0;
        }
      } end_f_loop(face, face_thread)
    }
    dynamic_thread = dynamic_thread->next;
  }

# endif /* !RP_HOST */
  Message("(%d) Synchronizing Fluent processes\n", process_id);
  PRF_GSYNC();
}

double get_ratio_edge_is_covered (
  face_t  face,
  Thread* face_thread,
  real    wave_height )
{
  real centroid_coords[ND_ND];
  real area_vector[ND_ND];
  real area = 0.0;
  real centroid_y = 0.0;
  real upper_y = 0.0;
  real lower_y = 0.0;

  F_CENTROID(centroid_coords, face, face_thread);
  F_AREA(area_vector, face, face_thread);
  area = NV_MAG(area_vector);
  centroid_y = centroid_coords[1];
  upper_y = centroid_y + area/2.0;
  lower_y = centroid_y - area/2.0;
  if (lower_y > surface_level + wave_height){
     return 0.0;
  }
  else if (upper_y < surface_level + wave_height){
     return 1.0;
  }
  return (surface_level + wave_height - lower_y) / area;
}

real get_ratio_triangle_is_covered (
  real coords[6],
  real wave_height )
{
  real ratio = 0.0;
  real min = coords[1];
  real max = coords[1];
  int above = 0;
  int below = 0;
  int above_count = 0;
  int i, j;
  real coords_leg_vertices[4];
  real lambda;
  real x_intersected_1;
  real x_intersected_2;
  real coords_small[6];
  real area, area_small;

  for (i=3; i < 6; i+=2){
    if (min > coords[i]) min = coords[i];
    if (max < coords[i]) max = coords[i];
  }

  if (min >= wave_height) return 0.0;
  if (max <= wave_height) return 1.0;

  /* Determine, whether one or two vertices are above interface */
  for (i=0; i < 3; i++){
    if (coords[i*2 +1] > wave_height){
      above_count++;
      above = i;
    }
    else {
      below = i;
    }
  }
  /*printf("above_count = %d, above = %d, below = %d\n", above_count, above, below);*/

  j=0;
  if (above_count == 1){ /* Select leg vertices below interface */
    coords_small[0] = coords[above*2];
    coords_small[1] = coords[above*2+1];
    for (i=0; i < 3; i++){
      if (i != above){
        coords_leg_vertices[j++] = coords[i*2];
        coords_leg_vertices[j++] = coords[i*2+1];
      }
    }
  }
  else { /* Select leg vertices above interface */
    coords_small[0] = coords[below*2];
    coords_small[1] = coords[below*2+1];
    for (i=0, j=0; i < 3; i++){
      if (i != below){
        coords_leg_vertices[j++] = coords[i*2];
        coords_leg_vertices[j++] = coords[i*2+1];
      }
    }
  }

  /*printf("coords head: %f, %f\n", coords_small[0], coords_small[1]);
  printf("leg vertex 1: %f, %f\n", coords_leg_vertices[0], coords_leg_vertices[1]);
  printf("leg vertex 2: %f, %f\n", coords_leg_vertices[2], coords_leg_vertices[3]);*/

  /* compute intersection points of interface with legs */
  lambda = (wave_height - coords_small[1]) / (coords_leg_vertices[1] - coords_small[1]);
  coords_small[2] = coords_small[0] + lambda*(coords_leg_vertices[0] - coords_small[0]);
  /*printf("lambda leg 1: %f, x leg 1: %f\n", lambda, coords_small[2]);*/
  lambda = (wave_height - coords_small[1]) / (coords_leg_vertices[3] - coords_small[1]);
  coords_small[4] = coords_small[0] + lambda*(coords_leg_vertices[2] - coords_small[0]);
  /*printf("lambda leg 2: %f, x leg 2: %f\n", lambda, coords_small[4]);*/
  coords_small[3] = wave_height;
  coords_small[5] = wave_height;

  area = compute_triangle_area_2d(coords);
  area_small = compute_triangle_area_2d(coords_small);
  /*printf("area = %f, area_small = %f\n", area, area_small);*/
  if (above_count == 1){
    return (area - area_small) / area;
  }
  return area_small / area;
}

real compute_triangle_area_2d(real coords[6])
{
  real ab[] = {coords[2] - coords[0], coords[3] - coords[1]};
  real ac[] = {coords[4] - coords[0], coords[5] - coords[1]};
  return fabs(0.5 * (ab[0]*ac[1] - ab[1]*ac[0]));
}
