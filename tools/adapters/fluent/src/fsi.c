#include "fsi.h"
#include "../../../../../src/precice/adapters/c/SolverInterfaceC.h"
#include "../../../../../src/precice/adapters/c/Constants.h"
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

double timestep_limit = 0.0;
double* forces = NULL;
int* force_indices = NULL;
int skip_grid_motion = BOOL_TRUE;
int did_gather_write_positions = BOOL_FALSE;
int did_gather_read_positions = BOOL_FALSE;
int thread_index = 0;
int dynamic_thread_size = 0;
int wet_edges_size = 0;
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
int comm_size = -1;
int require_create_checkpoint = BOOL_FALSE;
int* precice_force_ids; /* Gathered in host node (or serial node) */
int* precice_displ_ids;
/*int global_displ_index = 0;*/

#if ND_ND == 2
#define norm(a, b) sqrt((a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1]))
#else
#error Not implemented!
#endif

void count_dynamic_threads();
void gather_write_positions();
/*void gather_precice_write_indices();*/
void write_forces();
void gather_read_positions(Dynamic_Thread* dt);
void read_displacements(Dynamic_Thread* dt);
int check_write_positions();
int check_read_positions(Dynamic_Thread* dt);
void regather_read_positions(Dynamic_Thread* dt, int this_thread_size);
void regather_write_positions(int current_size);

void fsi_init(Domain* domain)
{
  int precice_process_id = -1; /* Process ID given to preCICE */
  printf("\nEntering fsi_init\n");

  #if !PARALLEL
  precice_process_id = 0;
  comm_size = 1;
  #else /* ! PARALELL*/
  #if RP_HOST
  precice_process_id = 0;
  #elif RP_NODE
  precice_process_id = myid + 1;
  #endif /* elif RP_NODE */
  comm_size = compute_node_count + 1;
  #endif /* else ! PARALLEL */

  Message("  (%d) Creating solver interface\n", myid);
  precicec_createSolverInterface("Fluent", precicec_nameConfiguration(),
                                precice_process_id, comm_size);

  Message("  (%d) Initializing coupled simulation\n", myid);
  timestep_limit = precicec_initialize();
  Message("  (%d) ... done\n", myid);

  #if !RP_HOST
  /* There might be several face threads forming the wet surface */
  count_dynamic_threads();
  #endif /* ! RP_HOST */

  if (precicec_isActionRequired(precicec_actionReadSimulationCheckpoint())){
    #if !RP_NODE /* HOST or SERIAL */
    Message("  (%d) Reading simulation checkpoint required\n", myid);
    RP_Set_Integer("udf/checkpoint", BOOL_TRUE);
    #endif /* !RP_NODE */
    did_gather_write_positions = BOOL_TRUE;
    did_gather_read_positions = BOOL_TRUE;
    skip_grid_motion = BOOL_FALSE; /* Read local displacements not stored in fluent checkpoint */
    precicec_fulfilledAction(precicec_actionReadSimulationCheckpoint());
  }

  if (precicec_isActionRequired(precicec_actionWriteIterationCheckpoint())){
    Message("  (%d) Implicit coupling\n", myid);
    #if ! RP_NODE
    RP_Set_Integer("udf/convergence", BOOL_FALSE);
    RP_Set_Integer("udf/iterate", BOOL_TRUE);
    #endif /* ! RP_NODE */
    precicec_fulfilledAction(precicec_actionWriteIterationCheckpoint());
  }
  else {
    Message("  (%d) Explicit coupling\n", myid);
  }

  Message("  (%d) Synchronizing Fluent processes\n", myid);
  PRF_GSYNC();
  printf("(%d) Leaving INIT\n", myid);
}

void fsi_write_and_advance()
{
  printf("(%d) Entering ON_DEMAND(write_and_andvance)\n", myid);
  int ongoing;
  int subcycling = ! precicec_isWriteDataRequired(CURRENT_TIMESTEP);
  int current_size = -1;

  /*Message("  (%d) write_and_advance 1\n", myid);*/
  if (subcycling){
    Message("  (%d) In subcycle, skip writing\n", myid);
  }
  else {
    if (!did_gather_write_positions){
      Message("  (%d) Gather write positions\n", myid);
      gather_write_positions();
    }
    else {
      /* Write positions can change in parallel mode when load-balancing occurs */
      current_size = check_write_positions();
      if (current_size != -1){
        regather_write_positions(current_size);
      }
    }
    #if !RP_HOST
    if (wet_edges_size){
      write_forces();
    }
    #endif
  }
  /*Message("  (%d) write_and_advance 2\n", myid);*/

  timestep_limit = precicec_advance(CURRENT_TIMESTEP);

  /* Read coupling state */
  ongoing = precicec_isCouplingOngoing();
  #if !RP_NODE
  RP_Set_Integer("udf/ongoing", ongoing);
  #endif /* !RP_NODE */

  if (precicec_isActionRequired(precicec_actionWriteIterationCheckpoint())){
    #if !RP_NODE
    RP_Set_Integer("udf/convergence", BOOL_TRUE);
    #endif /* !RP_NODE */
    precicec_fulfilledAction(precicec_actionWriteIterationCheckpoint());
  }

  if (precicec_isActionRequired(precicec_actionReadIterationCheckpoint())){
    #if !RP_NODE
    RP_Set_Integer("udf/convergence", BOOL_FALSE);
    #endif /* !RP_NODE */
    precicec_fulfilledAction(precicec_actionReadIterationCheckpoint());
  }

  #if !RP_NODE
  if (! precicec_isCouplingOngoing()){
    RP_Set_Integer("udf/convergence", BOOL_TRUE);
  }
  #endif /* !RP_NODE */

  if (precicec_isActionRequired(precicec_actionWriteSimulationCheckpoint())){
    precicec_fulfilledAction(precicec_actionWriteSimulationCheckpoint());
    require_create_checkpoint = BOOL_TRUE;
    #if !RP_NODE
    Message("  (%d) Writing simulation checkpoint required\n", myid);
    RP_Set_Integer("udf/checkpoint", BOOL_TRUE);
    #endif /* !RP_NODE */
  }
  #if !RP_NODE
  else {
    RP_Set_Integer("udf/checkpoint", BOOL_FALSE);
  }
  #endif /* !RP_NODE */

  printf("(%d) Leaving ON_DEMAND(write_and_advance)\n", myid);
}

void fsi_grid_motion(Domain* domain, Dynamic_Thread* dt, real time, real dtime)
{
  printf("\n(%d) Entering GRID_MOTION\n", myid);
  int meshID = precicec_getMeshID("WetSurface");
  int current_thread_size = -1;

  #if !RP_HOST /* Serial or node */
  if (thread_index == dynamic_thread_size){
    printf ("   Reset thread index\n");
    thread_index = 0;
  }
  printf("  (%d) Thread index = %d\n", myid, thread_index);
  Thread* face_thread  = DT_THREAD(dt);

  if (strncmp("gridmotions", dt->profile_udf_name, 11) != 0){
    printf("  (%d) ERROR: called gridmotions for invalid dynamic thread: %s\n",
            myid, dt->profile_udf_name);
    exit(1);
  }
  if (face_thread == NULL){
    printf("  (%d) ERROR: face_thread == NULL\n", myid);
    exit(1);
  }
  if (!did_gather_read_positions){
    gather_read_positions(dt);
  }
  else {
    /* Read positions can change in parallel mode when load-balancing occurs */
    current_thread_size = check_read_positions(dt);
    if (current_thread_size != -1){
      regather_read_positions(dt, current_thread_size);
    }
  }
  #endif /* !RP_HOST */

  if (skip_grid_motion){
    if (thread_index >= dynamic_thread_size-1){
      skip_grid_motion = BOOL_FALSE;
    }
    thread_index++;
    printf("  (%d) Skipping first round grid motion\n", myid);
    return;
  }

  /* Here the code was before */

  #if !RP_HOST
  SET_DEFORMING_THREAD_FLAG(THREAD_T0(face_thread));
  #endif /* !RP_HOST */

  precicec_mapReadData(meshID); /* Collective call necessary */

  #if !RP_HOST /* Serial or node */
  read_displacements(dt);
  thread_index++;
  #endif /* !RP_HOST */

  #if !RP_NODE
  Message("  (%d) convergence=%d, iterate=%d, couplingOngoing=%d\n",
          myid, RP_Get_Integer("udf/convergence"), RP_Get_Integer("udf/iterate"),
          precicec_isCouplingOngoing());
  if (RP_Get_Integer("udf/convergence") && RP_Get_Integer("udf/iterate") && precicec_isCouplingOngoing()){
    RP_Set_Integer("udf/convergence", BOOL_FALSE);
  }
  #endif /* !RP_NODE */
  if (! precicec_isCouplingOngoing()){
    precicec_finalize();
  }

  printf("(%d) Leaving GRID_MOTION\n", myid);
}

void fsi_plot_coords()
{
  printf("(%d) Entering ON_DEMAND(plot_coords)\n", myid);

  #if !RP_HOST
  int i=0, n=0;
  Domain* domain = NULL;
  Dynamic_Thread* dynamic_thread = NULL;
  Thread* face_thread = NULL;
  face_t face;
  Node* node = NULL;

  domain = Get_Domain(1);
  if (domain == NULL){
    printf("  (%d) ERROR: domain == NULL\n", myid);
    fflush(stdout);
    exit(1);
  }
  if (domain->dynamic_threads == NULL){
    printf("  (%d) ERROR: domain.dynamic_threads == NULL\n", myid);
    fflush(stdout);
    exit(1);
  }
  dynamic_thread = domain->dynamic_threads;
  while (dynamic_thread != NULL){
    if (strncmp("gridmotions", dynamic_thread->profile_udf_name, 11) == 0){
      face_thread  = DT_THREAD(dynamic_thread);
      if (face_thread == NULL){
        printf("  (%d) ERROR: face_thread == NULL\n", myid);
        fflush(stdout);
        exit(1);
      }
      begin_f_loop (face, face_thread){
        if (PRINCIPAL_FACE_P(face,face_thread)){
          f_node_loop (face, face_thread, n){
            node = F_NODE ( face, face_thread, n );
            if ((NODE_MARK(node) != 0) && ((NODE_MARK(node) != 1234))){
              printf("  (%d) ERROR: Unexpected node mark!\n", myid);
              fflush(stdout);
              exit(1);
            }
            if (NODE_MARK(node) != 1234){
              NODE_MARK(node) = 1234;
              /*if (i < 2){*/
                /*Message("coords: %.16E, %.16E\n", NODE_COORD(node)[0], NODE_COORD(node)[1]);*/
                /*fflush(stdout);*/
              /*}*/
              i++;
            }
          }
        }
      } end_f_loop(face, face_thread);
    }
    dynamic_thread = dynamic_thread->next;
  }

  /* Reset node mark */
  while (dynamic_thread != NULL){
    if (strncmp("gridmotions", dynamic_thread->profile_udf_name, 11) == 0){
      face_thread  = DT_THREAD(dynamic_thread);
      if (face_thread == NULL){
        printf("  (%d) ERROR: face_thread == NULL\n", myid);
        fflush(stdout);
        exit(1);
      }
      begin_f_loop (face, face_thread){
        if (PRINCIPAL_FACE_P(face,face_thread)){
          f_node_loop (face, face_thread, n){
            node = F_NODE ( face, face_thread, n );
            if (NODE_MARK(node) == 1234){
              NODE_MARK(node) = 0;
            }
          }
        }
      } end_f_loop(face, face_thread);
    }
    dynamic_thread = dynamic_thread->next;
  }
  #endif /* ! RP_HOST */

  printf("(%d) Leaving ON_DEMAND(plot_coords)\n", myid);
}

void count_dynamic_threads()
{
  printf("(%d) Entering count_dynamic_threads()\n", myid);
  Domain *domain = NULL;
  Dynamic_Thread* dynamic_thread = NULL;
  Thread* face_thread = NULL;
  face_t face;
  int node_index, i=0;
  Node* node = NULL;

  Message( "  (%d) counting dynamic threads: ", myid);
  domain = Get_Domain(1);
  if (domain == NULL){
    Message("  (%d) ERROR: domain == NULL\n", myid);
    exit(1);
  }
  if (domain->dynamic_threads == NULL){
    Message("  (%d) ERROR: domain.dynamic_threads == NULL\n", myid);
    exit (1);
  }
  dynamic_thread = domain->dynamic_threads;
  while (dynamic_thread != NULL){
    if (strncmp("gridmotions", dynamic_thread->profile_udf_name, 11) == 0){
      face_thread  = DT_THREAD(dynamic_thread);
      begin_f_loop (face, face_thread){ /* Thread face loop */
        if (PRINCIPAL_FACE_P(face,face_thread)){
          f_node_loop (face, face_thread, node_index){ /* Face node loop */
            node = F_NODE(face, face_thread, node_index);
            NODE_MARK(node) = 11111;
          }
        }
      } end_f_loop(face, face_thread)
      dynamic_thread_size++;
    }
    dynamic_thread = dynamic_thread->next;
  }
  dynamic_thread_node_size = (int*) malloc(dynamic_thread_size * sizeof(int));
  for (i=0;  i < dynamic_thread_size; i++){
    dynamic_thread_node_size[i] = 0;
  }

  /* Reset node marks */
  dynamic_thread = domain->dynamic_threads;
  while (dynamic_thread != NULL){
    if (strncmp("gridmotions", dynamic_thread->profile_udf_name, 11) == 0){
      face_thread  = DT_THREAD(dynamic_thread);
      begin_f_loop (face, face_thread){ /* Thread face loop */
        f_node_loop (face, face_thread, node_index){ /* Face node loop */
          node = F_NODE(face, face_thread, node_index);
          NODE_MARK(node) = 0;
        }
      } end_f_loop(face, face_thread)
    }
    dynamic_thread = dynamic_thread->next;
  }
  printf("  (%d) ... %d\n", myid, dynamic_thread_size);
  printf("(%d) Leaving count_dynamic_threads()\n", myid);
}

void gather_write_positions()
{
  printf("(%d) Entering gather_write_positions()\n", myid);
  #if !RP_HOST
  int meshID = precicec_getMeshID("WetSurface");
  int i = 0;
  double center[ND_ND];
  Domain* domain = NULL;
  Dynamic_Thread* dynamic_thread = NULL;
  Thread* face_thread = NULL;
  int thread_counter = 0;
  face_t face;

  Message("  (%d) Counting wet edges...\n", myid);
  domain = Get_Domain(1);
  if (domain == NULL){
    Message("  (%d) ERROR: domain == NULL\n", myid);
    exit(1);
  }
  if (domain->dynamic_threads == NULL){
    Message("  (%d) ERROR: domain.dynamic_threads == NULL\n", myid);
    exit(1);
  }
  dynamic_thread = domain->dynamic_threads;
  thread_counter = 0;
  while (dynamic_thread != NULL){
    if (strncmp("gridmotions", dynamic_thread->profile_udf_name, 11) == 0){
      Message("\n  (%d) Thread index %d\n", myid, thread_counter);
      face_thread  = DT_THREAD(dynamic_thread);
      if (face_thread == NULL){
        Message("  (%d) ERROR: face_thread == NULL\n", myid);
        exit(1);
      }
      begin_f_loop (face, face_thread){
        if (PRINCIPAL_FACE_P(face,face_thread)){
          wet_edges_size++;
        }
      } end_f_loop(face, face_thread);
      thread_counter++;
    }
    dynamic_thread = dynamic_thread->next;
  }
  Message("  (%d) ...done (counted %d wet edges)\n", myid, wet_edges_size);
  Message("  (%d) Allocating %d force vector values\n", myid, wet_edges_size * ND_ND);
  if (forces != NULL){
    Message("      (%d) ERROR: Forces vector allocated multiple times!\n", myid);
  }
  forces = (double*) malloc(wet_edges_size * ND_ND * sizeof(double));
  Message("  (%d) Allocating %d force indices\n", myid, wet_edges_size);
  force_indices = (int*) malloc(wet_edges_size * sizeof(int));
  /*force_indices_fluent = (int*) malloc(wet_edges_size * sizeof(int));*/
  Message("  (%d) Setting write positions...", myid);
  dynamic_thread = domain->dynamic_threads;
  thread_counter = 0;
  i = 0;
  while (dynamic_thread != NULL){
    if (strncmp("gridmotions", dynamic_thread->profile_udf_name, 11) == 0){
      Message("\n  (%d) Thread index %d\n", myid, thread_counter);
      face_thread  = DT_THREAD(dynamic_thread);
      if (face_thread == NULL){
        printf("  (%d) ERROR: face_thread == NULL\n", myid);
        fflush(stdout);
        exit(1);
      }
      begin_f_loop (face, face_thread){
        if (PRINCIPAL_FACE_P(face,face_thread)){
          F_CENTROID(center, face, face_thread);
          force_indices[i] = precicec_setWritePosition(meshID, center);
          F_UDMI(face, face_thread, 0) = force_indices[i];
          i++;
        }
      } end_f_loop(face, face_thread);
      thread_counter++;
    }
    dynamic_thread = dynamic_thread->next;
  }
  Message("  (%d) ...done counting wet edges\n", myid);
  did_gather_write_positions = BOOL_TRUE;
  #endif /* ! RP_HOST */

  /* Gather precice read and write indices */
  #if PARALLEL
  /*gather_precicec_write_indices();*/
  #endif /* PARALLEL */

  /* Setup precice index tables for checkpoint and load balancing */
  #if !RP_NODE /* Host or serial */
  #endif /* ! RP_NODE */
  printf("(%d) Leaving gather_write_positions()\n", myid);
}

void gather_read_positions(Dynamic_Thread* dt)
{
  printf("(%d) Entering gather_read_positions()\n", myid);
  Thread* face_thread  = DT_THREAD(dt);
  Node* node;
  face_t face;
  int n = 0, dim = 0;
  int array_index = 0, node_index = 0;
  double coords[ND_ND];
  int meshID = precicec_getMeshID("WetSurface");

  /* Count not yet as updated (from other threads) marked nodes */
  begin_f_loop(face, face_thread){
    if (PRINCIPAL_FACE_P(face,face_thread)){
      f_node_loop (face, face_thread, n){
        node = F_NODE ( face, face_thread, n );
        if (NODE_POS_NEED_UPDATE(node)){
          NODE_MARK(node) = 12345;
          wet_nodes_size++;
          dynamic_thread_node_size[thread_index]++;
        }
      }
    }
  } end_f_loop(face, face_thread);

  /* Get initial coordinates and reset update marking */
  printf("  (%d) Reallocating %d initial positions ...\n", myid, wet_nodes_size);
  initial_coords = (double*) realloc(initial_coords, wet_nodes_size * ND_ND * sizeof(double));
  displacements = (double*) realloc(displacements, wet_nodes_size * ND_ND * sizeof(double));
  displ_indices = (int*) realloc(displ_indices, wet_nodes_size * sizeof(int));
  array_index = wet_nodes_size - dynamic_thread_node_size[thread_index];
  begin_f_loop (face, face_thread){
    if (PRINCIPAL_FACE_P(face,face_thread)){
      f_node_loop(face, face_thread, n){
        node = F_NODE(face, face_thread, n);
        if (NODE_MARK(node) == 12345){
          NODE_MARK(node) = 1;  /*Set node to need update*/
          for (dim=0; dim < ND_ND; dim++){
            coords[dim] = NODE_COORD(node)[dim];
            initial_coords[array_index*ND_ND+dim] = coords[dim];
          }
          /*if (array_index < 10){
            printf("  (%d) initial coord %.16E\n", myid, initial_coords[array_index*ND_ND]);
            fflush(stdout);
          }*/
          node_index = precicec_setReadPosition(meshID, coords);
          displ_indices[array_index] = node_index;
          array_index++;
        }
      }
    }
  } end_f_loop(face, face_thread);
  printf("  (%d) Set %d (of %d) displacement read positions ...\n", myid,
          array_index - wet_nodes_size + dynamic_thread_node_size[thread_index],
          dynamic_thread_node_size[thread_index]);

  if (thread_index == dynamic_thread_size - 1){
    did_gather_read_positions = BOOL_TRUE;
  }
  printf("(%d) Leaving gather_read_positions()\n", myid);
}


void read_displacements(Dynamic_Thread* dt)
{
  int displID = precicec_getDataID("Displacements");
  int offset = 0;
  int i = 0, n = 0, dim = 0;
  Thread* face_thread  = DT_THREAD(dt);
  Node* node;
  face_t face;
  real max_displ_delta = 0.0;

  if (dynamic_thread_node_size[thread_index] > 0){
    Message("  (%d) Reading displacements...\n", myid);
    offset = 0;
    for (i = 0; i < thread_index; i++){
      offset += dynamic_thread_node_size[i];
    }
    precicec_readBlockVectorData(displID, dynamic_thread_node_size[thread_index],
        displ_indices + offset, displacements + ND_ND * offset);
    /* TEST TEST TEST */
    /*for (i=0; i < dynamic_thread_node_size[thread_index]; i++){
      displacements[ND_ND*offset + i*ND_ND] = 0.0;
      displacements[ND_ND*offset + i*ND_ND+1] = 0.0;
    }*/

    Message("  (%d) Setting displacements...\n", myid);
    i = offset * ND_ND;
    begin_f_loop (face, face_thread){
      if (PRINCIPAL_FACE_P(face,face_thread)){
        f_node_loop (face, face_thread, n){
          node = F_NODE(face, face_thread, n);
          if (NODE_POS_NEED_UPDATE(node)){
            NODE_POS_UPDATED(node);
            /*if (i < 4){
              printf("  (%d) init %.16E, %.16E\n", myid, initial_coords[i],
                     initial_coords[i+1]);
              printf("  (%d) displ %.16E, %.16E\n", myid, displacements[i],
                     displacements[i+1]);
              fflush(stdout);
            }*/
            for (dim=0; dim < ND_ND; dim++){
              NODE_COORD(node)[dim] = initial_coords[i+dim] + displacements[i+dim];
              /* displacements[i+dim] = NODE_COORD(node)[dim] - NODE_COORD_N(node)[dim];  store delta for rbf mesh motion */
              if (fabs(displacements[i+dim]) > fabs(max_displ_delta)){
                max_displ_delta = displacements[i + dim];
              }
            }
            i += ND_ND;
          }
        }
      }
    } end_f_loop (face, face_thread);
    Message("  (%d) ...done\n", myid);
  }
  Message("  (%d) Max displacement delta: %f\n", myid, max_displ_delta);
}

void write_forces()
{
  int forceID = precicec_getDataID("Forces");
  int i=0, j=0;
  Domain* domain = NULL;
  Dynamic_Thread* dynamic_thread = NULL;
  int thread_counter = 0;
  real area[ND_ND];
  real pressure_force[ND_ND];
  real viscous_force[ND_ND];
  double total_force[ND_ND];
  double max_force = 0.0;

  domain = Get_Domain(1);
  if (domain == NULL){
    Message("  (%d) ERROR: domain == NULL\n", myid);
    exit(1);
  }
  if (domain->dynamic_threads == NULL){
    Message("  (%d) ERROR: domain.dynamic_threads == NULL\n", myid);
    exit(1);
  }

  dynamic_thread = domain->dynamic_threads;
  thread_counter = 0;
  Message("  (%d) Gather forces...\n", myid);
  i = 0;
  while (dynamic_thread != NULL){
    if (strncmp("gridmotions", dynamic_thread->profile_udf_name, 11) == 0){
      Message("  (%d) Thread index %d\n", myid, thread_counter);
      Thread* face_thread  = DT_THREAD(dynamic_thread);
      if (face_thread == NULL){
        Message("  (%d) ERROR: face_thread == NULL\n", myid);
        exit(1);
      }
      face_t face;
      begin_f_loop (face, face_thread){
        if (PRINCIPAL_FACE_P(face,face_thread)){
          F_AREA(area, face, face_thread);
          NV_VS(viscous_force, =, F_STORAGE_R_N3V(face,face_thread,SV_WALL_SHEAR),*,-1.0);
          NV_VS(pressure_force, =, area, *, F_P(face,face_thread));
          NV_VV(total_force, =, viscous_force, +, pressure_force);
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
  Message("  (%d) ...done (with %d force values)\n", myid, i);
  Message("  (%d) Writing forces...\n", myid);
  precicec_writeBlockVectorData(forceID, wet_edges_size, force_indices, forces);
  Message("  (%d) ...done\n", myid );
  Message("  (%d) Max force: %f\n", max_force);
  if (thread_counter != dynamic_thread_size){
    Message ( "  (%d) ERROR: Number of dynamic threads has changed to %d!\n", myid, thread_counter );
    exit(1);
  }
}

int check_write_positions()
{
  #if !RP_HOST
  Domain* domain = NULL;
  Dynamic_Thread* dynamic_thread = NULL;
  Thread* face_thread = NULL;
  int thread_counter = 0;
  face_t face;
  int wet_edges_check_size = 0;

  Message("  (%d) Checking write positions...\n", myid);
  domain = Get_Domain(1);
  if (domain == NULL){
    Message("  (%d) ERROR: domain == NULL\n", myid);
    exit(1);
  }
  if (domain->dynamic_threads == NULL){
    Message("  (%d) ERROR: domain.dynamic_threads == NULL\n", myid);
    exit(1);
  }
  dynamic_thread = domain->dynamic_threads;
  thread_counter = 0;
  while (dynamic_thread != NULL){
    if (strncmp("gridmotions", dynamic_thread->profile_udf_name, 11) == 0){
      /*Message("\n  (%d) Thread index %d\n", myid, thread_counter);*/
      face_thread  = DT_THREAD(dynamic_thread);
      if (face_thread == NULL){
        Message("  (%d) ERROR: Thread %d: face_thread == NULL\n", myid, thread_counter);
        exit(1);
      }
      begin_f_loop (face, face_thread){
        if (PRINCIPAL_FACE_P(face,face_thread)){
          wet_edges_check_size++;
        }
      } end_f_loop(face, face_thread);
      thread_counter++;
    }
    dynamic_thread = dynamic_thread->next;
  }
  Message("  (%d) ...done (currently %d wet edges, old is %d)", myid,
          wet_edges_check_size, wet_edges_size);
  if (wet_edges_check_size != wet_edges_size) {
    return wet_edges_check_size;
  }
  #endif /* ! RP_HOST */
  return -1;
}

int check_read_positions(Dynamic_Thread* dt)
{
  Message("  (%d) Checking read positions...\n", myid);
  int i = 0, n = 0;
  Thread* face_thread = DT_THREAD(dt);
  Node* node;
  face_t face;

  /* Count nodes */
  begin_f_loop(face, face_thread){
    if (PRINCIPAL_FACE_P(face,face_thread)){
      f_node_loop (face, face_thread, n){
        node = F_NODE(face, face_thread, n);
        if (NODE_POS_NEED_UPDATE(node)){
          NODE_MARK(node) = 12345;
          i++;
        }
      }
    }
  } end_f_loop(face, face_thread);

  /* Reset node marks */
  begin_f_loop(face, face_thread){
    if (PRINCIPAL_FACE_P(face,face_thread)){
      f_node_loop (face, face_thread, n){
        node = F_NODE(face, face_thread, n);
        if (NODE_MARK(node) == 12345){
          NODE_MARK(node) = 1; /* Set node to need update*/
        }
      }
    }
  } end_f_loop(face, face_thread);

  if (i != dynamic_thread_node_size[thread_index]){
    Message("  (%d) Wet node count has changed for dynamic thread %d!\n",
            myid, thread_index);
    return i;
  }
  return -1;
}

void regather_read_positions(Dynamic_Thread* dt, int thread_new_size)
{
  Message("  (%d) Regathering read positions...\n", myid);

  int i=0, j=0, n=0, dim=0;
  Thread* face_thread = DT_THREAD(dt);
  Node* node;
  face_t face;
  int meshID = precicec_getMeshID("WetSurface");
  int all_size = precicec_getReadNodesSize(meshID);
  int displID = precicec_getDataID("Displacements");
  double* new_coords = (double*) malloc(thread_new_size * ND_ND * sizeof(double));
  int* new_indices = (int*) malloc(thread_new_size * sizeof(double));
  double* all_coords = (double*) malloc(all_size * ND_ND * sizeof(double));
  double* all_displ = (double*) malloc(all_size * ND_ND * sizeof(double));
  int* all_indices = (int*) malloc(all_size * sizeof(int));
  double* tail_coords = NULL;
  int* tail_indices = NULL;
  /*double coords[ND_ND];*/
  int left_i, right_i;
  int front_size, tail_size, new_size;

  for (i=0; i < all_size; i++){
    all_indices[i] = i;
  }
  precicec_readBlockVectorData(displID, all_size, all_indices, all_displ);
  precicec_getReadPositions(meshID, all_size, all_indices, all_coords);
  for (i=0; i < all_size*ND_ND; i++){
    if (i < all_size*ND_ND) {
      printf("  (%d) coods %.16E\n", myid, all_coords[i]);
      fflush(stdout);
    }
    all_coords[i] += all_displ[i];
  }

  /* Determine new indices, positions, and displacements */
  i = 0;
  begin_f_loop(face, face_thread){
    if (PRINCIPAL_FACE_P(face,face_thread)){
      f_node_loop (face, face_thread, n){
        node = F_NODE(face, face_thread, n);
        if (NODE_POS_NEED_UPDATE(node)){
          NODE_MARK(node) = 12345;
          for (j=0; j < all_size; j++){ /* Find position index in all positions */
            for (dim=0; dim < ND_ND; dim++){ /* Vector equality comp */
              if (i < 10){
                printf("  (%d) %.16E == %.16E\n", myid, NODE_COORD(node)[dim],
                       all_coords[j*ND_ND+dim]);
                fflush(stdout);
              }

              if (fabs(NODE_COORD(node)[dim] - all_coords[j*ND_ND+dim]) > 1e-7){
                break;
              }
            }
            if (dim == ND_ND){ /* If equal */
              /*printf("  (%d) Equal!\n", myid); fflush(stdout);*/
              new_indices[i] = all_indices[j];
              for (dim=0; dim < ND_ND; dim++){
                left_i = i*ND_ND+dim;
                right_i = j*ND_ND+dim;
                new_coords[left_i] = all_coords[right_i] - all_displ[right_i];
              }
              break;
            }
          }
          if (j == all_size){
            printf("  (%d) ERROR: Didn't find suitable index while regathering read indices!\n", myid);
            fflush(stdout);
            exit(1);
          }
          i++;
        }
      }
    }
  } end_f_loop(face, face_thread);

  /* Reset node marks */
  begin_f_loop(face, face_thread){
    if (PRINCIPAL_FACE_P(face,face_thread)){
      f_node_loop (face, face_thread, n){
        node = F_NODE(face, face_thread, n);
        if (NODE_MARK(node) == 12345){
          NODE_MARK(node) = 1; /* Set node to need update*/
        }
      }
    }
  } end_f_loop(face, face_thread);

  /* Count entries before this face thread and after it */
  front_size = 0;
  for (i=0; i < thread_index; i++){
    front_size += dynamic_thread_node_size[i];
  }
  tail_size = 0;
  for (i=thread_index+1; i < dynamic_thread_size; i++){
    tail_size += dynamic_thread_node_size[i];
  }
  tail_coords = (double*) malloc(tail_size*ND_ND*sizeof(double));
  tail_indices = (int*) malloc(tail_size*sizeof(int));

  /* Store tail entries */
  for (i=0; i < tail_size; i++){
    j = front_size + dynamic_thread_node_size[thread_index] + i;
    for (dim=0; dim < ND_ND; dim++){
      tail_coords[i*ND_ND+dim] = initial_coords[j*ND_ND+dim];
    }
    tail_indices[i] = displ_indices[j];
  }

  /* Insert new and tail entries */
  new_size = front_size + thread_new_size + tail_size;
  initial_coords = (double*) realloc(initial_coords, new_size * ND_ND * sizeof(double));
  displacements = (double*) realloc(displacements, new_size * ND_ND * sizeof(double));
  displ_indices = (int*) realloc(displ_indices, new_size * sizeof(int));
  for (i=0; i < thread_new_size; i++){
    j = front_size + i;
    for (dim=0; dim < ND_ND; dim++){
      initial_coords[j*ND_ND+dim] = new_coords[i*ND_ND+dim];
    }
    printf("  (%d) new index: %d\n", myid, new_indices[i]);
    displ_indices[j] = new_indices[i];
  }
  for (i=0; i < tail_size; i++){
    j = front_size + thread_new_size + i;
    for (dim=0; dim < ND_ND; dim++){
      initial_coords[j*ND_ND+dim] = tail_coords[i*ND_ND+dim];
    }
    displ_indices[j] = tail_indices[i];
  }

  wet_nodes_size -= dynamic_thread_node_size[thread_index];
  wet_nodes_size += thread_new_size;
  dynamic_thread_node_size[thread_index] = thread_new_size;

  free(all_indices);
  free(all_coords);
  free(new_indices);
  free(new_coords);
  free(tail_indices);
  free(tail_coords);

  Message("  (%d) ... done regathering read positions...\n", myid); fflush(stdout);
}

void regather_write_positions(int current_size)
{
  #if !RP_HOST
  int i = 0;
  Domain* domain = NULL;
  Dynamic_Thread* dynamic_thread = NULL;
  Thread* face_thread = NULL;
  int thread_counter = 0;
  face_t face;

  Message("  (%d) Regather write positions...\n", myid);

  forces = (double*) realloc(forces, current_size * ND_ND * sizeof(double));
  force_indices = (int*) realloc(force_indices, current_size * sizeof(int));

  domain = Get_Domain(1);
  if (domain == NULL){
    Message("  (%d) ERROR: domain == NULL\n", myid);
    exit(1);
  }
  if (domain->dynamic_threads == NULL){
    Message("  (%d) ERROR: domain.dynamic_threads == NULL\n", myid);
    exit(1);
  }
  dynamic_thread = domain->dynamic_threads;
  thread_counter = 0;
  i = 0;
  while (dynamic_thread != NULL){
    if (strncmp("gridmotions", dynamic_thread->profile_udf_name, 11) == 0){
      Message("\n  (%d) Thread index %d\n", myid, thread_counter);
      face_thread  = DT_THREAD(dynamic_thread);
      if (face_thread == NULL){
        Message("  (%d) ERROR: face_thread == NULL\n", myid);
        exit(1);
      }
      begin_f_loop (face, face_thread){
        if (PRINCIPAL_FACE_P(face,face_thread)){
          force_indices[i] = F_UDMI(face, face_thread, 0);
          i++;
        }
      } end_f_loop(face, face_thread);
      thread_counter++;
    }
    dynamic_thread = dynamic_thread->next;
  }
  wet_edges_size = current_size;
  Message("  (%d) ...done regathering write positions\n", myid);
  #endif /* ! RP_HOST */
}
