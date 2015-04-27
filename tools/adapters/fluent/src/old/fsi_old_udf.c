#include "udf.h"
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

double timestep_limit = 0.0;
char* const createIterationCheckpoint = "write-iteration-checkpoint";
char* const readIterationCheckpoint = "read-iteration-checkpoint";
double* forces = NULL;
int* force_indices = NULL;
int* force_indices_fluent = NULL;
int skip_grid_motion = BOOL_TRUE;
int gather_write_positions = BOOL_TRUE;
int gather_read_positions = BOOL_TRUE;
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
int precice_process_id = -1;
int comm_size = -1;
int require_create_checkpoint = BOOL_FALSE;
/*int global_displ_index = 0;*/

#if ND_ND == 2
#define norm(a, b) sqrt((a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1]))
#else
#error Not implemented!
#endif

void lu_decomposition(double* matrix, int* pivots, int n);

/**
 * Assumes ones at the diagonal of matrix.
 */
void forward_substitution(const double* matrix, const double* b, double* x, int n);
void backward_substitution(const double* matrix, const double* b, double* x, int n);
void export_matrix(const double* matrix, int n, const char* name);
void export_vector(const double* vector, int n, const char* name);
void export_int_vector(const int* vector, int n, const char* name);
void compute_rbf_matrix();
void compute_rbf_coefficients();

DEFINE_INIT(init,d)
{
  Dynamic_Thread* dynamic_thread = NULL;
  Thread* face_thread = NULL;
  face_t face;
  int node_index;
  Node* node = NULL;
  int i = 0; /* Loop counter */
  int dim = 0; /* Dimension loop counter */
  Message("\nEntering INIT\n");

# if RP_NODE
  precice_process_id = myid + 1;
  comm_size = compute_node_count + 1;
  Message("   Compute node: id=%d, size=%d\n", myid, comm_size);
# endif
# if RP_HOST
  precice_process_id = 0;
  comm_size = compute_node_count + 1;
  Message("   Host: id=%d, size=%d\n", myid, comm_size);
# endif
# if !PARALLEL
  precice_process_id = 0;
  comm_size = 1;
  Message("   (%d) Serial node: id=%d, size=%d\n", myid, precice_process_id, comm_size);
# endif

  Message("   (%d) Creating coupling interface\n", myid);
  precice_createSolverInterface("Fluent", precice_nameConfiguration(),
                                precice_process_id, comm_size);

  Message("   (%d) Initializing coupled simulation\n", myid);
  timestep_limit = precice_initialize();
  Message("   (%d) ... done\n", myid);

  if (precice_isActionRequired(precice_actionReadSimulationCheckpoint())){
#   if !RP_NODE
    Message("   (%d) Reading simulation checkpoint required\n", myid);
    RP_Set_Integer("udf/checkpoint", BOOL_TRUE);
#   endif /* !RP_NODE */
    gather_write_positions = BOOL_FALSE;
    gather_read_positions = BOOL_FALSE;
    /*precice_fulfilledAction(precice_actionReadSimulationCheckpoint());*/
  }

  if (precice_isActionRequired(createIterationCheckpoint)){
    Message("   (%d) Implicit coupling\n", myid);
#   if !RP_NODE
    RP_Set_Integer("udf/convergence", BOOL_FALSE);
    RP_Set_Integer("udf/iterate", BOOL_TRUE);
#   endif /* !RP_NODE */
    precice_fulfilledAction(createIterationCheckpoint);
  }
  else {
    Message("   (%d) Explicit coupling\n", myid);
  }

# if !RP_HOST
  /* Count dynamic threads */
  Message ( "   (%d) counting dynamic threads: ", myid );
  Domain *domain = Get_Domain(1);
  if (domain == NULL){
    Message("      (%d) ERROR: domain == NULL\n", myid);
    exit(1);
  }
  if (domain->dynamic_threads == NULL){
    Message ( "      (%d) ERROR: domain.dynamic_threads == NULL\n", myid );
    exit (1);
  }
  dynamic_thread = domain->dynamic_threads;
  while (dynamic_thread != NULL){
    if (strncmp("gridmotions", dynamic_thread->profile_udf_name, 11) == 0){
      face_thread  = DT_THREAD(dynamic_thread);
      begin_f_loop (face, face_thread){ /* Thread face loop */
        f_node_loop (face, face_thread, node_index){ /* Face node loop */
          node = F_NODE(face, face_thread, node_index);
          /*Message("   (%d) dynamic thread node mark = %d\n", myid, NODE_MARK(node));*/
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

  /* Comment in when doing own mesh movement
  Message ( "   (%d) counting boundary nodes: ", myid );
  thread_loop_f (face_thread, domain){ / Domain thread loop /
    if (THREAD_T1(face_thread) == NULL){
      begin_f_loop (face, face_thread){ / Thread face loop /
        f_node_loop (face, face_thread, node_index){ / Face node loop /
          node = F_NODE(face, face_thread, node_index);
          /Message("   (%d) node mark = %d\n", myid, NODE_MARK(node));/
          if ((NODE_MARK(node) != 12345) && (NODE_MARK(node) != 11111)){
            NODE_MARK(node) = 12345;
            boundary_nodes_size++;
          }
        }
      } end_f_loop(face, face_thread)
    }
  }
  Message("%d\n", boundary_nodes_size);
  */

  /*
  Message ( "   (%d) counting deformable nodes: ", myid );
  thread_loop_f (face_thread, domain){
    if (THREAD_T1(face_thread) != NULL){
      begin_f_loop (face, face_thread){
        f_node_loop (face, face_thread, node_index){
          node = F_NODE(face, face_thread, node_index);
          if (NODE_MARK(node) != 12346){
            NODE_MARK(node) = 12346;
            deformable_nodes_size++;
          }
        }
      } end_f_loop(face, face_thread)
    }
  }
  Message("%d\n", deformable_nodes_size);
  */

  /*
  boundary_coords = (double*) malloc(boundary_nodes_size * ND_ND * sizeof(double));
  Message("   (%d) Storing boundary node coords...\n", myid);
  thread_loop_f (face_thread, domain){ / Domain thread loop /
    if (THREAD_T1(face_thread) == NULL){
      begin_f_loop (face, face_thread){ / Thread face loop /
        f_node_loop (face, face_thread, node_index){ / Face node loop /
          node = F_NODE(face, face_thread, node_index);
          if (NODE_MARK(node) == 12345){
            NODE_MARK(node) = 0;  / ??? not sure, should reset node /
            for (dim=0; dim < ND_ND; dim++){
              boundary_coords[i+dim] = NODE_COORD(node)[dim];
            }
            /Message("      (%d) Boundary coords[%d]: %f, %f\n", myid, i,
                boundary_coords[i], boundary_coords[i+1]);/
            i += ND_ND;
          }
        }
      } end_f_loop(face, face_thread)
    }
  }
  Message("      (%d) ...done", myid);
  */

  dynamic_thread = domain->dynamic_threads;
  while (dynamic_thread != NULL){
    if (strncmp("gridmotions", dynamic_thread->profile_udf_name, 11) == 0){
      face_thread  = DT_THREAD(dynamic_thread);
      begin_f_loop (face, face_thread){ /* Thread face loop */
        f_node_loop (face, face_thread, node_index){ /* Face node loop */
          node = F_NODE(face, face_thread, node_index);
          /*Message("   (%d) reset dynamic thread node mark = %d\n", myid, NODE_MARK(node));*/
          NODE_MARK(node) = 0;
        }
      } end_f_loop(face, face_thread)
    }
    dynamic_thread = dynamic_thread->next;
  }

# endif /* !RP_HOST */
  Message("(%d) Synchronizing Fluent processes\n", myid);
  PRF_GSYNC();
  Message("(%d) Leaving INIT\n", myid);
}

DEFINE_ON_DEMAND(write_and_advance)
{
  Message("(%d) Entering ON_DEMAND(write_and_andvance)\n", myid);
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
      Message("      (%d) ERROR: domain == NULL\n", myid);
      exit(1);
    }
    /*Message ( "      Domain name = %s", domain->name );*/
    if (domain->dynamic_threads == NULL){
      Message("      (%d) ERROR: domain.dynamic_threads == NULL\n", myid);
      exit(1);
    }

    if (gather_write_positions){
      Message("      (%d) Counting wet edges...", myid);
      dynamic_thread = domain->dynamic_threads;
      thread_counter = 0;
      while (dynamic_thread != NULL){
        if (strncmp("gridmotions", dynamic_thread->profile_udf_name, 11) == 0){
          Message("\n      (%d) Thread index %d\n", myid, thread_counter);
          Thread* face_thread  = DT_THREAD(dynamic_thread);
          if (face_thread == NULL){
            Message("      (%d) ERROR: face_thread == NULL\n", myid);
            exit(1);
          }
          face_t face;
          begin_f_loop (face, face_thread){
            if (PRINCIPAL_FACE_P(face,face_thread)){
              wet_edges_size++;
            }
          } end_f_loop(face, face_thread);
          thread_counter++;
        }
        dynamic_thread = dynamic_thread->next;
      }
      Message("      (%d) ...done (counted %d wet edges)", myid, wet_edges_size);
      Message("      (%d) Allocating %d force vector values\n", myid, wet_edges_size * ND_ND);
      if (forces != NULL){
        Message("      (%d) ERROR: Forces vector allocated multiple times!\n", myid);
      }
      forces = (double*) malloc(wet_edges_size * ND_ND * sizeof(double));
      Message("      (%d) Allocating %d force indices\n", myid, wet_edges_size);
      force_indices = (int*) malloc(wet_edges_size * sizeof(int));
      force_indices_fluent = (int*) malloc(wet_edges_size * sizeof(int));

      Message("      (%d) Setting write positions...", myid);
      dynamic_thread = domain->dynamic_threads;
      thread_counter = 0;
      i = 0;
      while (dynamic_thread != NULL){
        if (strncmp("gridmotions", dynamic_thread->profile_udf_name, 11) == 0){
          Message("\n      (%d) Thread index %d\n", myid, thread_counter);
          Thread* face_thread  = DT_THREAD(dynamic_thread);
          if (face_thread == NULL){
            Message("      (%d) ERROR: face_thread == NULL\n", myid);
            exit(1);
          }
          face_t face;
          begin_f_loop (face, face_thread){
            if (PRINCIPAL_FACE_P(face,face_thread)){
              F_CENTROID(center, face, face_thread);
              force_indices[i] = precice_setWritePosition(meshID, center);
              /*cell = F_C0(face,face_thread);*/
              /*force_indices_fluent[i] =*/
              i++;
            }
          } end_f_loop(face, face_thread);
          thread_counter++;
        }
        dynamic_thread = dynamic_thread->next;
      }
      Message("      (%d) ...done\n", myid);
      gather_write_positions = BOOL_FALSE;
    }
    else {
      /* Check if thread has changed */
      Message("\n      (%d) Checking dynamic threads for changes\n", myid);
      dynamic_thread = domain->dynamic_threads;
      thread_counter = 0;
      i = 0;
      while (dynamic_thread != NULL){
        if (strncmp("gridmotions", dynamic_thread->profile_udf_name, 11) == 0){
          Message("\n      (%d) Thread index %d\n", myid, thread_counter);
          Thread* face_thread  = DT_THREAD(dynamic_thread);
          if (face_thread == NULL){
            Message("      (%d) ERROR: face_thread == NULL\n", myid);
            exit(1);
          }
          face_t face;
          begin_f_loop (face, face_thread){
            if (PRINCIPAL_FACE_P(face,face_thread)){

              i++;
            }
          } end_f_loop(face, face_thread);
          thread_counter++;
        }
        dynamic_thread = dynamic_thread->next;
      }
      Message("      (%d) ...done\n", myid);
      if (i != wet_edges_size){
        Message("\n      (%d) Wet edge count has changed to %d\n", myid, i);
        /*
        face = C_FACE(cell,cell_thread,local_face_index);
        cell = F_C0(face,face_thread);
        */
        /*
        Message("\n      (%d) Re-indexing force locations\n", myid);
        force_indices[i] = precice_setWritePosition(meshID, center);
        Message("      (%d) ...done\n", myid);
        */
      }
    }

    dynamic_thread = domain->dynamic_threads;
    thread_counter = 0;
    Message("      (%d) Gather forces...\n", myid);
    i = 0;
    while (dynamic_thread != NULL){
      if (strncmp("gridmotions", dynamic_thread->profile_udf_name, 11) == 0){
        Message("\n      (%d) Thread index %d\n", myid, thread_counter);
        Thread* face_thread  = DT_THREAD(dynamic_thread);
        if (face_thread == NULL){
          Message("      (%d) ERROR: face_thread == NULL\n", myid);
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
            /*Message ( "      (%d) Area = %f, %f; Visc.Force = %f, %f; PressureForce = %f, %f; Total = %f, %f\n",
                      myid, area[0], area[1], viscous_force[0], viscous_force[1],
                      pressure_force[0], pressure_force[1], total_force[0], total_force[1] );
            */
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
    Message("      (%d) ...done (with %d force values)\n", myid, i);
    Message("      (%d) Writing forces...\n", myid);
    precice_writeBlockVectorData(forceID, wet_edges_size, force_indices, forces);
    Message("      (%d) ...done\n", myid );
    Message("      (%d) Max force: %f\n", max_force);
    if (thread_counter != dynamic_thread_size){
      Message ( "      (%d) ERROR: Number of dynamic threads has changed to %d!\n", myid, thread_counter );
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
    Message("   (%d) Writing simulation checkpoint required\n", myid);
    RP_Set_Integer("udf/checkpoint", BOOL_TRUE);
#   endif /* !RP_NODE */
  }
# if !RP_NODE
  else {
    RP_Set_Integer("udf/checkpoint", BOOL_FALSE);
  }
# endif /* !RP_NODE */


  Message("(%d) Leaving ON_DEMAND(write_and_advance)\n", myid);
}

DEFINE_GRID_MOTION(gridmotions,domain,dt,time,dtime)
{
  Message("\n(%d) Entering GRID_MOTION\n", myid);
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
    Message("      (%d) Reallocating %d initial positions ...\n", myid, wet_nodes_size);
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
    Message("      (%d) Set %d (of %d) displacement read positions ...\n", myid,
            array_index - wet_nodes_size + dynamic_thread_node_size[thread_index],
            dynamic_thread_node_size[thread_index]);

    if (thread_index == dynamic_thread_size - 1){
      gather_read_positions = BOOL_FALSE;
      /*compute_rbf_matrix();*/
    }
  }
  else {
    i = 0;
    begin_f_loop(face, face_thread){
      if (PRINCIPAL_FACE_P(face,face_thread)){
        f_node_loop (face, face_thread, n){
          node = F_NODE(face, face_thread, n);
          if (NODE_POS_NEED_UPDATE(node)){
            i++;
          }
        }
      }
    } end_f_loop(face, face_thread);
    if (i != dynamic_thread_node_size[thread_index]){
      Message("\n      (%d) Wet node count has changed to %d\n", myid, i);
    }
  }

  /*indexNode = 0;*/
  SET_DEFORMING_THREAD_FLAG(THREAD_T0(face_thread));
# endif /* !RP_HOST */

  precice_mapReadData(meshID); /* Collective call necessary */

# if !RP_HOST /* Serial or node */
  if (dynamic_thread_node_size[thread_index] > 0){
    Message("      (%d) Reading displacements...\n", myid);
    offset = 0;
    for (i = 0; i < thread_index; i++){
      offset += dynamic_thread_node_size[i];
    }
    /*Message("      offset=%d, size=%d\n", offset, dynamic_thread_node_size[thread_index]);*/
    precice_readBlockVectorData(displID, dynamic_thread_node_size[thread_index],
        displ_indices + offset, displacements + ND_ND * offset);

    Message("      (%d) Setting displacements...\n", myid);
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
            /*Message("old coord = %f, %f", NODE_COORD_N(node)[0], NODE_COORD_N(node)[1] );
            Message("; displ = %f, %f", displacements[i], displacements[i+1]);
            */
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
            /*Message("; new coord = %f, %f\n", NODE_COORD(node)[0], NODE_COORD(node)[1]);*/
            i += ND_ND;
          }
        }
      }
    } end_f_loop (face, face_thread);
    Message("      (%d) ...done\n", myid);

    /*
    if (thread_index == dynamic_thread_size-1){
      Message("      (%d) Moved nodes counter: %d, of %d\n", myid,
              moved_nodes_counter, );
      compute_rbf_coefficients();
      moved_nodes_counter = 0;
    }
    */
  }
  Message("      (%d) Max displacement delta: %f\n", myid, max_displ_delta);

  thread_index++;
# endif /* !RP_HOST */

# if !RP_NODE
  Message("      (%d) convergence=%d, iterate=%d, couplingOngoing=%d\n",
          myid, RP_Get_Integer("udf/convergence"), RP_Get_Integer("udf/iterate"),
          precice_isCouplingOngoing());
  if (RP_Get_Integer("udf/convergence") && RP_Get_Integer("udf/iterate") && precice_isCouplingOngoing()){
    RP_Set_Integer("udf/convergence", BOOL_FALSE);
  }
# endif /* !RP_NODE */
  if (! precice_isCouplingOngoing()){
    precice_finalize();
  }

  Message("(%d) Leaving GRID_MOTION\n", myid);
}

DEFINE_GEOM(deformmesh, domain, dt, position)
{
  double* boundPos = NULL; /* For Fluent boundary coordinates */
  int i=0; /* Loop counter */
  double displacement[2]; /* Displacement value */
  double r = 0.0; /* Radial distance of position to boundary node */
  double coeff = 0.0; /* Radial basis function coefficient */

  if (initial_coords != NULL){
    displacement[0] = 0.0;
    displacement[1] = 0.0;
    for (i=0; i < wet_nodes_size + boundary_nodes_size; i++){
      boundPos = i < wet_nodes_size ? &(initial_coords[i*ND_ND])
                                    : &(boundary_coords[(i-wet_nodes_size)*ND_ND]);
      r = norm(position, boundPos);
      /* Thin plate spline radial basis function evaluation */
      coeff = r > 1e-10 ? log10(r) * pow(r, 2) : 0.0;
      displacement[0] += coeff * x_coeff_vector[i];
      displacement[1] += coeff * y_coeff_vector[i];
    }

    Message("(%d) Moving mesh node %f, %f by %f, %f\n", myid, position[0],
            position[1], displacement[0], displacement[1]);
    position[0] += displacement[0];
    position[1] += displacement[1];
    moved_nodes_counter++;
  }
}

DEFINE_RW_FILE(readdata, filep)
{
  int i = 0;
  int dim = 0;
  int compute_node_count_check = -1;

# if PARALLEL
  int rank = -1;
  int wet_nodes_size_buffer = -1;
  int wet_edges_size_buffer = -1;
  real* initial_coords_buffer = NULL;
  int* displ_indices_buffer = NULL;
  int* force_indices_buffer = NULL;
# endif /* PARALLEL */

  Message("\n (%d) Read data\n", myid);
  if (precice_isActionRequired(precice_actionReadSimulationCheckpoint())){
    Message("(%d) Reading data from file\n", myid);
#   if !PARALLEL
    /* Number of compute nodes */
    fscanf(filep, "%d ", &compute_node_count_check);
    if (compute_node_count != compute_node_count_check){
      Message("(%d) ERROR: Compute node count inconsistent process number (checkpt: %d, now: %d)!\n",
              myid, compute_node_count_check, compute_node_count);
      exit(1);
    }
    /* Number of nodes holding displacements */
    fscanf(filep, "%d", &wet_nodes_size);
    Message("(%d) Number of wet nodes: %d\n", myid, wet_nodes_size);
    dynamic_thread_node_size[0] = wet_nodes_size; /* HACK */
    Message("(%d) Allocating memory for nodes ...", myid);
    initial_coords = (double*) malloc(wet_nodes_size * ND_ND * sizeof(double));
    displacements = (double*) malloc(wet_nodes_size * ND_ND * sizeof(double));
    displ_indices = (int*) malloc(wet_nodes_size * sizeof(int));
    Message(" done.\n");
    /* Initial coordinates displacement nodes */
    Message("(%d) Read initial coords: ", myid);
    for (i=0; i < wet_nodes_size; i++){
      fscanf(filep, "%lf", initial_coords + i*ND_ND);
      Message("%f", initial_coords[i*ND_ND]);
      for (dim=1; dim < ND_ND; dim++){
        fscanf(filep, "%lf", initial_coords + i*ND_ND + dim);
        Message(", %f", initial_coords[i*ND_ND + dim]);
      }
      Message("\n");
    }
    /* Displacement indices */
    Message("(%d) Read displ indices: ", myid);
    for (i=0; i < wet_nodes_size; i++){
      fscanf(filep, "%d", displ_indices + i);
      Message("%d ", displ_indices[i]);
    }
    Message("\n");
    /* Number of nodes holding forces */
    fscanf(filep, "%d", &wet_edges_size);
    Message("(%d) Number of wet edges: %d\n", myid, wet_edges_size);
    Message("(%d) Allocating memory for edges ...", myid);
    forces = (double*) malloc(wet_edges_size * ND_ND * sizeof(double));
    force_indices = (int*) malloc(wet_edges_size * sizeof(int));
    Message(" done.\n");
    /* Force indices */
    Message("(%d) Read force indices: ", myid);
    for (i=0; i < wet_edges_size; i++){
      fscanf(filep, "%d", force_indices + i);
      Message("%d ", force_indices[i]);
    }
    Message("\n");

#   else /* !PARALLEL -> PARALLEL */

#   if RP_HOST
    Message("(%d) Read as host\n", myid);
    /* Number of compute nodes */
    fscanf(filep, "%d ", &compute_node_count_check);
    if (compute_node_count != compute_node_count_check){
      Message("(%d) ERROR: Compute node count inconsistent (checkpoint: %d, current: %d)!\n",
              myid, compute_node_count_check, compute_node_count);
      exit(1);
    }
    compute_node_loop(rank){
      Message("(%d) Treating rank %d\n", myid, rank);
      /* Number of nodes holding displacements */
      fscanf(filep, "%d", &wet_nodes_size_buffer);
      Message("(%d) Number of wet nodes: %d\n", myid, wet_nodes_size_buffer);
      PRF_CSEND_INT(node_zero, &wet_nodes_size_buffer, 1, myid);

      if (wet_nodes_size_buffer > 0){
        initial_coords_buffer = (double*) malloc(wet_nodes_size_buffer * ND_ND * sizeof(double));
        displ_indices_buffer = (int*) malloc(wet_nodes_size_buffer * sizeof(int));

        Message("(%d) Read initial coords: ", myid);
        for (i=0; i < wet_nodes_size_buffer; i++){
          fscanf(filep, "%lf", initial_coords_buffer + i*ND_ND);
          Message("%f", initial_coords_buffer[i*ND_ND]);
          for (dim=1; dim < ND_ND; dim++){
            fscanf(filep, "%lf", initial_coords_buffer + i*ND_ND + dim);
            Message(", %f", initial_coords_buffer[i*ND_ND + dim]);
          }
          Message("; ");
        }
        PRF_CSEND_REAL(node_zero, initial_coords_buffer, wet_nodes_size_buffer*ND_ND, myid);

        Message("(%d) Read displ indices: ", myid);
        for (i=0; i < wet_nodes_size_buffer; i++){
          fscanf(filep, "%d", displ_indices_buffer + i);
          Message("%d ", displ_indices_buffer[i]);
        }
        Message("\n");
        PRF_CSEND_INT(node_zero, displ_indices_buffer, wet_nodes_size_buffer, myid);

        free(initial_coords_buffer);
        free(displ_indices_buffer);
      }

      fscanf(filep, "%d", &wet_edges_size_buffer);
      Message("(%d) Number of wet edges: %d\n", myid, wet_edges_size_buffer);
      PRF_CSEND_INT(node_zero, &wet_edges_size_buffer, 1, myid);

      if (wet_edges_size_buffer > 0){
        force_indices_buffer = (int*) malloc(wet_edges_size_buffer * sizeof(int));
        Message("(%d) Read force indices: ", myid);
        for (i=0; i < wet_edges_size_buffer; i++){
          fscanf(filep, "%d", force_indices_buffer + i);
          Message("%d ", force_indices_buffer[i]);
        }
        Message("\n");
        PRF_CSEND_INT(node_zero, force_indices_buffer, wet_edges_size_buffer, myid);
        free(force_indices_buffer);
      }
    }
#   endif /* RP_HOST */

#   if RP_NODE
    if (myid == node_zero){
      Message("\n (%d) Read as node zero\n", myid);
      PRF_CRECV_INT(node_host, &wet_nodes_size, 1, node_host);
      initial_coords = (double*) malloc(wet_nodes_size * ND_ND * sizeof(double));
      displacements = (double*) malloc(wet_nodes_size * ND_ND * sizeof(double));
      displ_indices = (int*) malloc(wet_nodes_size * sizeof(int));
      if (wet_nodes_size > 0){
        PRF_CRECV_REAL(node_host, initial_coords, wet_nodes_size*ND_ND, node_host);
        PRF_CRECV_INT(node_host, displ_indices, wet_nodes_size, node_host);
      }
      PRF_CRECV_INT(node_host, &wet_edges_size, 1, node_host);
      forces = (double*) malloc(wet_edges_size * ND_ND * sizeof(double));
      force_indices = (int*) malloc(wet_edges_size * sizeof(int));
      if (wet_edges_size > 0){
        PRF_CRECV_INT(node_host, force_indices, wet_edges_size, node_host);
      }
      compute_node_loop_not_zero(rank){
        PRF_CRECV_INT(node_host, &wet_nodes_size_buffer, 1, node_host);
        PRF_CSEND_INT(rank, &wet_nodes_size_buffer, 1, myid);
        if (wet_nodes_size_buffer > 0){
          initial_coords_buffer = (double*) malloc(wet_nodes_size_buffer * ND_ND * sizeof(double));
          displ_indices_buffer = (int*) malloc(wet_nodes_size_buffer * sizeof(int));
          PRF_CRECV_REAL(node_host, initial_coords_buffer, wet_nodes_size_buffer*ND_ND, node_host);
          PRF_CSEND_REAL(rank, initial_coords_buffer, wet_nodes_size_buffer*ND_ND, myid);
          PRF_CRECV_INT(node_host, displ_indices_buffer, wet_nodes_size_buffer, node_host);
          PRF_CSEND_INT(rank, displ_indices_buffer, wet_nodes_size_buffer, myid);
          free(initial_coords_buffer);
          free(displ_indices_buffer);
        }
        PRF_CRECV_INT(node_host, &wet_edges_size_buffer, 1, node_host);
        PRF_CSEND_INT(rank, &wet_edges_size_buffer, 1, myid);
        if (wet_edges_size_buffer > 0){
          force_indices_buffer = (int*) malloc(wet_edges_size_buffer * sizeof(int));
          PRF_CRECV_INT(node_host, force_indices_buffer, wet_edges_size_buffer, node_host);
          PRF_CSEND_INT(rank, force_indices_buffer, wet_edges_size_buffer, myid);
          free(force_indices_buffer);
        }
      }
    }
    else {
      Message("\n (%d) Read as other node\n", myid);
      PRF_CRECV_INT(node_zero, &wet_nodes_size, 1, node_zero);
      initial_coords = (double*) malloc(wet_nodes_size * ND_ND * sizeof(double));
      displacements = (double*) malloc(wet_nodes_size * ND_ND * sizeof(double));
      displ_indices = (int*) malloc(wet_nodes_size * sizeof(int));
      if (wet_nodes_size > 0){
        PRF_CRECV_REAL(node_zero, initial_coords, wet_nodes_size*ND_ND, node_zero);
        PRF_CRECV_INT(node_zero, displ_indices, wet_nodes_size, node_zero);
      }
      PRF_CRECV_INT(node_zero, &wet_edges_size, 1, node_zero);
      forces = (double*) malloc(wet_edges_size * ND_ND * sizeof(double));
      force_indices = (int*) malloc(wet_edges_size * sizeof(int));
      if (wet_edges_size > 0){
        PRF_CRECV_INT(node_zero, force_indices, wet_edges_size, node_zero);
        Message("(%d) Received force indices: ", myid);
        for (i=0; i < wet_edges_size; i++){
          Message("%d ", force_indices[i]);
        }
        Message("\n");
      }
    }
    dynamic_thread_node_size[0] = wet_nodes_size; /* HACK */
#   endif /* RP_NODE */

#   endif /* !PARALLEL */
    precice_fulfilledAction(precice_actionReadSimulationCheckpoint());
    Message("(%d) Done reading data\n", myid);
  }
}

DEFINE_RW_FILE(writedata, filep)
{
  int i = 0;
  int dim = 0;
  int appendcheckpoint = BOOL_FALSE;

# if PARALLEL
  int rank = -1;
  int wet_nodes_size_buffer = -1;
  int wet_edges_size_buffer = -1;
  real* initial_coords_buffer = NULL;
  int* displ_indices_buffer = NULL;
  int* force_indices_buffer = NULL;
# endif /* PARALLEL */

  Message("\n (%d) Write data\n", myid);
# if !RP_NODE
  appendcheckpoint = RP_Get_Integer("udf/appendcheckpoint");
# endif
  host_to_node_int_1(appendcheckpoint);

  if (require_create_checkpoint && appendcheckpoint){
    require_create_checkpoint = BOOL_FALSE;

#   if !PARALLEL
    Message("(%d) Writing data to file\n", myid);
    /* Number of compute nodes */
    fprintf(filep, "%d ", 1);
    /* Number of nodes holding displacements */
    fprintf(filep, "%d ", wet_nodes_size);
    /* Initial coordinates */
    Message("(%d) Writing initial coords:\n", myid);
    for (i=0; i < wet_nodes_size; i++){
      fprintf(filep, "%f", initial_coords[i*ND_ND]);
      Message("%f", initial_coords[i*ND_ND]);
      for (dim=1; dim < ND_ND; dim++){
        fprintf(filep, " %f", initial_coords[i*ND_ND + dim]);
        Message(" %f", initial_coords[i*ND_ND]);
      }
      fprintf(filep, " ");
      Message(" ");
    }
    /* Displacement indices */
    Message("(%d) Writing displ indices: ", myid);
    for (i=0; i < wet_nodes_size; i++){
      fprintf(filep, "%d ", displ_indices[i]);
      Message(" %f", displ_indices[i]);
    }
    Message("(%d) ... done.\n", myid);
    /* Number of nodes holding forces */
    fprintf(filep, "%d ", wet_edges_size);
    /* Force indices */
    for (i=0; i < wet_edges_size; i++){
      fprintf(filep, "%d ", force_indices[i]);
    }

#   else  /* !PARALLEL -> PARALLEL*/

#   if RP_HOST /* Gather infos from rank 0 and write to file */
    Message("\n (%d) Write as host\n", myid);
    fprintf(filep, "%d ", compute_node_count);
    compute_node_loop(rank){
      Message("Treating rank %d\n", rank);

      PRF_CRECV_INT(node_zero, &wet_nodes_size_buffer, 1, node_zero);
      Message("  Wet node size: %d\n", wet_nodes_size_buffer);
      fprintf(filep, "%d ", wet_nodes_size_buffer);

      if (wet_nodes_size_buffer > 0){
        /* Write initial coordinates */
        Message("  Initial coords: ");
        initial_coords_buffer = (real*) malloc(wet_nodes_size_buffer * ND_ND * sizeof(real));
        PRF_CRECV_REAL(node_zero, initial_coords_buffer, wet_nodes_size_buffer * ND_ND, node_zero);
        for (i=0; i < wet_nodes_size_buffer; i++){
          fprintf(filep, "%f", initial_coords_buffer[i*ND_ND]);
          Message("%f", initial_coords_buffer[i*ND_ND]);
          for (dim=1; dim < ND_ND; dim++){
            fprintf(filep, " %f", initial_coords_buffer[i*ND_ND + dim]);
            Message(" %f", initial_coords_buffer[i*ND_ND]);
          }
          fprintf(filep, " ");
          Message(" ");
        }
        Message("\n");
        free(initial_coords_buffer);

        /* Write discplacement indices */
        displ_indices_buffer = (int*) malloc(wet_nodes_size_buffer * sizeof(int));
        PRF_CRECV_INT(node_zero, displ_indices_buffer, wet_nodes_size_buffer, node_zero);
        Message("  Displacement indices:");
        for (i=0; i < wet_nodes_size_buffer; i++){
          fprintf(filep, "%d ", displ_indices_buffer[i]);
          Message(" %d", displ_indices_buffer[i]);
        }
        Message("\n");
        free(displ_indices_buffer);
      }

      /* Number of nodes holding forces */
      PRF_CRECV_INT(node_zero, &wet_edges_size_buffer, 1, node_zero);
      Message("  Wet edge size: %d\n", wet_edges_size_buffer);
      fprintf(filep, "%d ", wet_edges_size_buffer);

      if (wet_edges_size_buffer > 0){
        /* Force indices */
        force_indices_buffer = (int*) malloc(wet_edges_size_buffer * sizeof(int));
        PRF_CRECV_INT(node_zero, force_indices_buffer, wet_edges_size_buffer, node_zero);
        Message("  Force indices: ");
        for (i=0; i < wet_edges_size_buffer; i++){
          Message(" %d", force_indices_buffer[i]);
          fprintf(filep, "%d ", force_indices_buffer[i]);
        }
        free(force_indices_buffer);
      }
      Message("\n\n");
    }
#   endif /* RP_HOST */

#   if RP_NODE
    if (myid == node_zero){
      Message("\n (%d) Write as node zero\n", myid);
      /* Send own data to host */
      PRF_CSEND_INT(node_host, &wet_nodes_size, 1, myid);
      if (wet_nodes_size > 0){
        PRF_CSEND_REAL(node_host, initial_coords, wet_nodes_size * ND_ND, myid);
        PRF_CSEND_INT(node_host, displ_indices, wet_nodes_size, myid);
      }
      PRF_CSEND_INT(node_host, &wet_edges_size, 1, myid);
      if (wet_edges_size > 0){
        PRF_CSEND_INT(node_host, force_indices, wet_edges_size, myid);
      }

      /* Receive data of other compute nodes and forward it to host */
      compute_node_loop_not_zero(rank){
        PRF_CRECV_INT(rank, &wet_nodes_size_buffer, 1, rank);
        PRF_CSEND_INT(node_host, &wet_nodes_size_buffer, 1, myid);

        if (wet_nodes_size_buffer > 0){
          initial_coords_buffer = (real*) malloc(wet_nodes_size_buffer * ND_ND * sizeof(real));
          PRF_CRECV_REAL(rank, initial_coords_buffer, wet_nodes_size_buffer * ND_ND, rank);
          PRF_CSEND_REAL(node_host, initial_coords_buffer, wet_nodes_size_buffer * ND_ND, myid);
          free(initial_coords_buffer);

          displ_indices_buffer = (int*) malloc(wet_nodes_size_buffer * sizeof(int));
          PRF_CRECV_INT(rank, displ_indices_buffer, wet_nodes_size_buffer, rank);
          PRF_CSEND_INT(node_host, displ_indices_buffer, wet_nodes_size_buffer, myid);
          free(displ_indices_buffer);
        }

        PRF_CRECV_INT(rank, &wet_edges_size_buffer, 1, rank);
        PRF_CSEND_INT(node_host, &wet_edges_size_buffer, 1, myid);

        if (wet_edges_size_buffer > 0){
          force_indices_buffer = (int*) malloc(wet_edges_size_buffer * sizeof(int));
          PRF_CRECV_INT(rank, force_indices_buffer, wet_edges_size_buffer, rank);
          PRF_CSEND_INT(node_host, force_indices_buffer, wet_edges_size_buffer, myid);
          free(force_indices_buffer);
        }
      }
    }
    else {
      Message("\n (%d) Write as other node\n", myid);
      /* Send own data to node zero */
      PRF_CSEND_INT(node_zero, &wet_nodes_size, 1, myid);
      if (wet_nodes_size > 0){
        PRF_CSEND_REAL(node_zero, initial_coords, wet_nodes_size * ND_ND, myid);
        PRF_CSEND_INT(node_zero, displ_indices, wet_nodes_size, myid);
      }
      PRF_CSEND_INT(node_zero, &wet_edges_size, 1, myid);
      if (wet_edges_size > 0){
        PRF_CSEND_INT(node_zero, force_indices, wet_edges_size, myid);
      }
    }
#   endif /* RP_NODE */
#   endif /* !PARALLEL */
  }
}

void lu_decomposition(double* matrix, int* pivots, int n)
{
  int k=0, i=0, j=0;
  int maxIndex=0;
  double max=0, current=0, temp=0;
  for (k=0; k < n; k++){
    /* Determine line with max pivot */
    maxIndex = k;
    max = fabs(matrix[k*n + k]);
    for (i=k+1; i < n; i++){
      current = fabs(matrix[i*n + k]);
      if (current > max){
        maxIndex = i;
        max = current;
      }
    }
    pivots[k] = maxIndex;

    /* Exchange lines */
    for (j=0; j < n; j++){
      temp = matrix[k*n + j];
      matrix[k*n + j] = matrix[maxIndex*n + j];
      matrix[maxIndex*n + j] = temp;
    }

    /* Compute scaling elements */
    for (i=k+1; i < n; i++){
      matrix[i*n + k] /= matrix[k*n + k];
    }

    /* Subtract contributions from each line */
    for (i=k+1; i < n; i++){
      for (j=k+1; j < n; j++){
        matrix[i*n + j] -= matrix[i*n + k] * matrix[k*n + j];
      }
    }
  }
}

/**
 * Assumes ones at the diagonal of matrix.
 */
void forward_substitution(const double* matrix, const double* b, double* x, int n)
{
  int i=0, j=0;
  x[0] = b[0];
  for (i=1; i < n; i++){
    x[i] = b[i];
    for (j=0; j < i; j++){
      x[i] -= matrix[i*n + j] * x[j];
    }
  }
}

void backward_substitution(const double* matrix, const double* b, double* x, int n)
{
  int i=0, j=0;
  for (i=n-1; i >= 0; i--){
    x[i] = b[i];
    for (j=i+1; j < n; j++){
      x[i] -= matrix[i*n + j] * x[j];
    }
    x[i] /= matrix[i*n + i];
  }
}

void export_matrix(const double* matrix, int n, const char* name)
{
  int i=0, j=0;
  FILE* file;
  file = fopen(name, "w");
  for (i=0; i < n; i++){
    for (j=0; j < n; j++){
      fprintf(file, "%.16f ", matrix[i*n + j]);
    }
    fprintf(file, "\n");
  }
  fclose(file);
}

void export_vector(const double* vector, int n, const char* name)
{
  int i=0;
  FILE* file;
  file = fopen(name, "w");
  for (i=0; i < n; i++){
    fprintf(file, "%.16f\n", vector[i]);
  }
  fclose(file);
}

void export_int_vector(const int* vector, int n, const char* name)
{
  int i=0;
  FILE* file;
  file = fopen(name, "w");
  for (i=0; i < n; i++){
    fprintf(file, "%d\n", vector[i]);
  }
  fclose(file);
}

void compute_rbf_matrix()
{
  int entries_per_line = 0;
  int i=0, j=0, dim=0;
  double* x_i = NULL;
  double* x_j = NULL;
  double r = 0.0;
  Message("      (%d) Compute rbf mapping matrix...\n", myid);
  entries_per_line = wet_nodes_size + boundary_nodes_size + ND_ND + 1;
  Message("      (%d) ... allocate arrays (size=%d)\n", myid, entries_per_line);
  if (c_matrix != NULL){
    Message("      (%d) ERROR: compute_rbf_matrix called multiple times!\n", myid);
    exit(1);
  }
  c_matrix = (double*) malloc(entries_per_line * entries_per_line * sizeof(double));
  b_vector = (double*) malloc(entries_per_line * sizeof(double));
  x_coeff_vector = (double*) malloc(entries_per_line * sizeof(double));
  y_coeff_vector = (double*) malloc(entries_per_line * sizeof(double));
  pivots_vector = (int*) malloc(entries_per_line * sizeof(int));

  /* Compute C matrix */
  Message("      (%d) ... compute c matrix\n", myid);
  for (i=0; i < wet_nodes_size + boundary_nodes_size; i++){
    x_i = i < wet_nodes_size ? &(initial_coords[i*ND_ND]) : &(boundary_coords[(i-wet_nodes_size)*ND_ND]);
    for (j=i; j < wet_nodes_size + boundary_nodes_size; j++){
      x_j = j < wet_nodes_size ? &(initial_coords[j*ND_ND]) : &(boundary_coords[(j-wet_nodes_size)*ND_ND]);
      r = norm(x_i, x_j);
      /* Thin plate spline radial basis function evaluation */
      c_matrix[i*entries_per_line + j] = r > 1e-10 ? log10(r) * pow(r, 2) : 0.0;
    }
    c_matrix[(i+1)*entries_per_line - ND_ND - 1] = 1.0;
    for (dim=0; dim < ND_ND; dim++){
      c_matrix[(i+1)*entries_per_line - ND_ND + dim] = x_i[dim];
    }
  }
  /* Fill lower right corner with zeros */
  for (i=0; i < ND_ND+1; i++){
    for (j=i; j < ND_ND+1; j++){
      c_matrix[(entries_per_line-ND_ND-1+i)*entries_per_line + wet_nodes_size + boundary_nodes_size + j] = 0.0;
    }
  }
  /* Copy symmetric lower diagonal part */
  Message("      (%d) ... mirror symmetric part of c matrix\n", myid);
  for (i=0; i < entries_per_line; i++){
    for (j=0; j < entries_per_line; j++){
      c_matrix[j*entries_per_line + i] = c_matrix[i*entries_per_line + j];
    }
  }
  Message("      (%d) ... export c matrix\n", myid);
  export_matrix(c_matrix, entries_per_line, "c_matrix.txt");
  Message("      (%d) ... lu-decompose c matrix\n", myid);
  lu_decomposition(c_matrix, pivots_vector, entries_per_line);
  export_matrix(c_matrix, entries_per_line, "c_lu_matrix.txt");
  export_int_vector(pivots_vector, entries_per_line, "pivots_vector.txt");
  Message("      (%d) ... done\n", myid);
}

void compute_rbf_coefficients()
{
  Message("      (%d) Computing mesh deformation coefficients...\n", myid);
  int n = wet_nodes_size + boundary_nodes_size + ND_ND + 1;
  double temp = 0.0;
  int i = 0, j = 0;
  double* temp_vector = (double*) malloc(n * sizeof(double));
  /* Assign righthand-side vector values for x-displacement */
  for (i=0; i < wet_nodes_size; i++){
    b_vector[i] = displacements[i*2];
  }
  for (i=0; i < boundary_nodes_size; i++){
    b_vector[wet_nodes_size+i] = 0.0;
  }
  for (i=0; i < ND_ND+1; i++){
    b_vector[wet_nodes_size + boundary_nodes_size + i] = 0.0;
  }
  export_vector(b_vector, n, "vector_b_x_displ.txt");

  /* Compute coefficients for x-displacement */
  for (i=0; i < n; i++){
    temp = b_vector[i];
    b_vector[i] = b_vector[pivots_vector[i]];
    b_vector[pivots_vector[i]] = temp;
  }
  forward_substitution(c_matrix, b_vector, temp_vector, n);
  backward_substitution(c_matrix, temp_vector, x_coeff_vector, n);
  export_vector(x_coeff_vector, n, "vector_coeff_x_displ.txt");

  /* Assign righthand-side vector values for x-displacement */
  for (i=0; i < wet_nodes_size; i++){
    b_vector[i] = displacements[i*2 + 1];
  }
  for (i=0; i < boundary_nodes_size; i++){
    b_vector[wet_nodes_size + i + 1] = 0.0;
  }
  for (i=0; i < ND_ND+1; i++){
    b_vector[wet_nodes_size + boundary_nodes_size + i] = 0.0;
  }
  export_vector(b_vector, n, "vector_b_y_displ.txt");

  /* Compute coefficients for y-displacement */
  for (i=0; i < n; i++){
    temp = b_vector[i];
    b_vector[i] = b_vector[pivots_vector[i]];
    b_vector[pivots_vector[i]] = temp;
  }
  forward_substitution(c_matrix, b_vector, temp_vector, n);
  backward_substitution(c_matrix, temp_vector, y_coeff_vector, n);
  export_vector(y_coeff_vector, n, "vector_coeff_y_displ.txt");

  free(temp_vector);
  Message("      (%d) ... done\n", myid);
}
