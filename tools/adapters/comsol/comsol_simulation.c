#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
// ---- these includes are related to the
#include <glib.h>
#include <rpc/rpc.h>
#include <sys/types.h>
#include <netinet/in.h>
#include <netdb.h>

#include "scriptext.h"
#include "mpi.h"
#include "fsi_mesh.h"
#include "fsi_map-mat.h"
#include "fsi_interface_socket.h"

#define dist(x1,x2,y1,y2)         sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2))
#define abs(x)                    (x>0)?x:-x
#define epsilon                   0.000000000000001

#define SIM_LOAD              "Benchmark_init"
#define SIM_FORCES            "Benchmark_new_forces"
#define SIM_NEXT_STEP         "Benchmark_new_step"
#define SIM_NEXT_STEP_RESTART "Benchmark_new_step_restart"
#define SIM_RESTART           "Benchmark_restart"

#define COMMUNICATION_PORT 53219
#define HOSTNAME           "localhost"
#define SOCKET_ERROR       -1
#define QUEUE_SIZE         5

/**
 * Debug print macros.
 *
 * If CPP flag "Debug" is set, debug output is written into files.
 */
#ifdef Debug

#define comsolDebug(filepointer, ...) \
        fprintf(filepointer, __VA_ARGS__);

#define comsolOpen(filepointer, filename, accesstype) \
        filepointer = fopen(filename, accesstype);

#define comsolClose(filepointer) \
        fclose(filepointer);

#else

#define comsolDebug(filepointer, ...)
#define comsolOpen(filepointer, filename, accesstype)
#define comsolClose(filepointer)

#endif

static int mesh_point_count = 0;
static int x_mesh_point_count = 0;
static int y_mesh_point_count = 0;
static int tmp_mesh_point_count = 0;
static int point_count = 0;
static int init_points = 0;
static int iteration_number=0;
static int *coresp_x = NULL;
static int *coresp_y = NULL;
// this is the mesh
static FSI_Mesh* mesh = NULL;
static FSI_Data* mesh_forces = NULL;
static FSI_Data* mesh_displacements = NULL;
static FSI_Data* mesh_displacementdeltas = NULL;
static FSI_Data* mesh_velocities = NULL;
static FSI_Data* mesh_velocitydeltas = NULL;
//static FILE *f = NULL;
static FILE* f_debug = NULL;

/** This variable should be recieved from the supervisor */
static double supervisor_timestep = 0.01;
/** This is the communication socket, Comsol waits always waits for the other side !!!*/
static int socketToPrec;
static int CommSocket;
/** This integer will show if we make n new step or we "remake" the step */
static int NewStep = 1;
static int repeatStep = 1;

// The function for setting the X forces
CL_EXPORT void force_x
(
  clEnv*  env,
  int     nOut,
  clData* out[],
  int     nIn,
  clData* in[] )
{
  //printf("ENTERING force_x\n"); fflush(stdout);
  int i = 0;
  int j;
  int ind;
  size_t size_elems_x = 0;
  size_t size_elems_y = 0;
  double * x = NULL;
  double * y = NULL;
  double * val = NULL;
  double * out_val = NULL;
  double min_dist;
  double tmp_real;

  //comsolOpen(f_debug, "comsol_debug.tmp", "a+");
  // get the coordinates
  size_elems_x = clGetNElems(in[0]);
  size_elems_y = clGetNElems(in[1]);
  x = clGetRealPtr(in[0]);
  y = clGetRealPtr(in[1]);

  //comsolDebug(f_debug, " force_x: size_elems_x = %d, size_elems_y = %d, x_mesh_point_count = %d\n",
  //            size_elems_x, size_elems_y, x_mesh_point_count );
  // set the output
  out[0] = clNewFull2D(env, CL_REAL, 1, size_elems_x);
  out_val = clGetRealPtr(out[0]);
	if (mesh_point_count == 0){
		//for (i=0; i < size_elems_x; i++){
    //  out_val[i] = 0;
    //}
		// this is only for this case when we want to count the points
		tmp_mesh_point_count += size_elems_x;
	}
	//else {
  for (i=0; i < size_elems_x; i++){
    if (coresp_x[x_mesh_point_count] >= 0){
      out_val[i] = FSI_Data_get_value(mesh_forces, coresp_x[x_mesh_point_count++], 1);
    }
    else {
      // get the most apropiate point for this point
      min_dist = dist(x[i],FSI_Mesh_get_node_x(mesh,0),y[i],FSI_Mesh_get_node_y(mesh,0));
      ind = 0;
      for (j = 1; j < FSI_Mesh_get_num_nodes(mesh); j++){
        tmp_real = dist(x[i] , FSI_Mesh_get_node_x(mesh,j) , y[i] , FSI_Mesh_get_node_y(mesh,j));
        if (tmp_real < min_dist){
          min_dist = tmp_real;
          ind = j;
        }
      }
      coresp_x[x_mesh_point_count] = ind;
      out_val[i] = FSI_Data_get_value(mesh_forces, coresp_x[x_mesh_point_count++], 1);
    }
    //comsolDebug(f_debug, "x-force %d = %f\n", i, out_val[i]);
  }
	//}
	//printf("LEAVING force_x\n"); fflush(stdout);
  return;
}

// The function for setting the Y forces
CL_EXPORT void force_y
(
  clEnv*  env,
  int     nOut,
  clData* out[],
  int     nIn,
  clData* in[] )
{
  //printf("ENTERING force_y\n"); fflush(stdout);
  int i = 0 , j , ind ;
  size_t size_elems_x = 0, size_elems_y = 0;
  double * x = NULL, * y = NULL;
  double * val = NULL, * out_val = NULL, min_dist , tmp_real;
  //comsolOpen(f_debug, "comsol_debug.tmp","a+");
  // get the coordinates
  size_elems_x = clGetNElems(in[0]);
  size_elems_y = clGetNElems(in[1]);
  x = clGetRealPtr(in[0]);
  y = clGetRealPtr(in[1]);
  //comsolDebug(f_debug," force_y %d %d : %d\n", size_elems_x ,size_elems_y , y_mesh_point_count);
  // set the output
  out[0] = clNewFull2D (env, CL_REAL, 1, size_elems_x);
  out_val = clGetRealPtr (out[0]);
	if (mesh_point_count == 0){
		//for (i=0; i<size_elems_x; i++){
		//  out_val[i] = 0;
		//}
		// this is only for this case when we want to count the points
		tmp_mesh_point_count += size_elems_x;
	}
	//else {
  for (i=0; i<size_elems_x; i++){
    if ( coresp_y[y_mesh_point_count] >= 0 ){
            //fprintf(f_debug," force_y found %d %d \n" , y_mesh_point_count , coresp_y[y_mesh_point_count]);
      /*fprintf(f_debug," Y Coord : %f %f : %d : %f %f\n" ,
          FSI_Mesh_get_node_x(mesh,coresp_y[y_mesh_point_count]) ,
          FSI_Mesh_get_node_y(mesh,coresp_y[y_mesh_point_count]) ,
          coresp_y[y_mesh_point_count] ,
          x[i] , y[i] );
      */
      out_val[i] = FSI_Data_get_value( mesh_forces , coresp_y[y_mesh_point_count++] , 2 );
      //fprintf(f_debug," Y Force set %f \n", out_val[i] );
    }
    else {
      // get the most apropiate point for this point
      min_dist = dist(x[i],FSI_Mesh_get_node_x(mesh,0),y[i],FSI_Mesh_get_node_y(mesh,0));
      ind = 0;
      for ( j = 1 ; j < FSI_Mesh_get_num_nodes(mesh) ; j++){
        tmp_real = dist(x[i] , FSI_Mesh_get_node_x(mesh,j) , y[i] , FSI_Mesh_get_node_y(mesh,j));
        if ( tmp_real < min_dist){
          min_dist = tmp_real;
          ind = j;
        }
      }

      //fprintf(f_debug," force_y %d %d  - %f %f Dist:%f - %f,%f \n" ,
      //     y_mesh_point_count , ind,  x[i],y[i], min_dist ,
      //     FSI_Mesh_get_node_x(mesh,ind) ,
      //     FSI_Mesh_get_node_y(mesh,ind) );
      coresp_y[y_mesh_point_count] = ind;
      out_val[i] = FSI_Data_get_value( mesh_forces ,
          coresp_y[y_mesh_point_count++] , 2 );
      //fprintf(f_debug," Y Force set %f \n", out_val[i] );
    }
    //comsolDebug(f_debug, "y-force %d = %f\n", i, out_val[i]);
  }

	//}
  //printf("LEAVING force_y\n"); fflush(stdout);
  return;
}

// this function is called after a new iteration
// it is needed to reset the point counters
CL_EXPORT void nextstep_simulation()
//(
//  clEnv *  env,
//  int      nOut,
//  clData * out[],
//  int      nIn,
//  clData * in[] )
{
	int i;
	//comsolOpen(f_debug, "comsol_debug.tmp","a+");
	comsolDebug(f_debug," nextstep called ... %d , %d \n" , mesh_point_count , tmp_mesh_point_count);

	x_mesh_point_count = 0;
	y_mesh_point_count = 0;
	// we should create this structure only at the beginning of the simualtion
	if (mesh_point_count == 0){
		mesh_point_count = tmp_mesh_point_count;
		comsolDebug(f_debug," nextstep called ... create corespond %d \n", mesh_point_count);
		//if (coresp_x != NULL){
		//    g_free(coresp_x);
		//    g_free(coresp_y);
		//}
		if (coresp_x == NULL){
		  coresp_x = g_malloc(100000*sizeof(int));
		  coresp_y = g_malloc(100000*sizeof(int));
		}
		for ( i = 0 ; i <   100000 ; i++){
			coresp_x[i] = -1;
			coresp_y[i] = -1;
		}
	}
	comsolDebug(f_debug," END nextstep called ... \n");
	//comsolClose(f_debug);
	tmp_mesh_point_count = 0;
	return;
}
// the cleanup function, to delete allocated memory
CL_EXPORT void cleanup_simulation (clEnv *env, int nOut, clData *out[], int nIn, clData *in[]) {
	//comsolOpen(f_debug, "comsol_debug.tmp","a+");
	comsolDebug(f_debug," ceanup called ... \n");
	//comsolClose(f_debug);

	FSI_Mesh_Data_destroy( mesh, "Forces" );
	FSI_Mesh_Data_destroy( mesh, "Displacements" );
	FSI_Mesh_Data_destroy( mesh, "Velocities" );
	FSI_Mesh_Data_destroy( mesh, "VelocityDeltas" );
	FSI_Mesh_destroy( mesh );
	mesh = NULL;
	mesh_forces = NULL;
	mesh_point_count = 0;
	x_mesh_point_count = 0;
	y_mesh_point_count = 0;
	tmp_mesh_point_count = 0;
	point_count = 0;
	init_points = 0;
	g_free(coresp_x);
	g_free(coresp_y);
	coresp_x = NULL;
	coresp_y = NULL;
	return;
}

// this function returns the iteration number which currently is not used
CL_EXPORT void return_itnumber (clEnv *env, int nOut, clData *out[], int nIn, clData *in[]) {
   out[0] = clNewReal (env, (double)iteration_number);
	//comsolOpen(f_debug,"comsol_debug.tmp","a+");
	comsolDebug(f_debug,"Return Itnumber  ... %d \n",iteration_number);
	//comsolClose(f_debug);
   return;
}

// the function which returns to the COMSOL side the imposed time step by the supervisor
CL_EXPORT void return_timestep (clEnv *env, int nOut, clData *out[], int nIn, clData *in[]) {
   out[0] = clNewReal (env, (double)supervisor_timestep);
	//comsolOpen(f_debug,"comsol_debug.tmp","a+");
	comsolDebug(f_debug,"Return Timestep  ... %f \n", supervisor_timestep );
	//comsolClose(f_debug);
	return;
}

// this function tells if we need to change the order of the points
int revertDirection(double x1, double y1 , double x2 , double y2 , double x3 , double y3){
  double m;
  double tmp,b;

  if ( (x2 - x1)*(x2 - x1) < epsilon){
     if ( (x1 - x3)*(y2-y1) > epsilon )
       return 0;// means do not revert
     else
       return 1;
  }else
  if ((y2 - y1)*(y2 - y1) < epsilon){
     if ( (y1 - y3)*(x2-x1) > epsilon )
       return 1;
     else
       return 0;// means do not revert
  }
  else{ // this is a simple line
     m = (y2 - y1)/(x2 - x1);
     b = y1 - m*x1;
     tmp = m*x3 + b - y3;
     if ( tmp * (x2-x1) < epsilon )
        return 0;// means do not revert
     else
        return 1;
  }
  return 0;
}

/** This function saves the actual workspace (state of the structure simulation)*/
void save_checkpoint(clEnv *env,const char *name){
  printf("Writing checkpoint!\n");
  comsolDebug(f_debug," Saving checkpoint ...\n", name);
  //const char command[50] = "save ";
  //strcpy(&(command[5]),name);
  // just evaluate the command
  // this command should evaluate in the workspace the following line:
  // "save checkpointfilename" then a "checkpointfilename.ws" file will be saved
  clEvalExpr(env,"flsave comsol_checkpoint fem",0);
  //clEvalExpr(env,"save comsol_checkpoint.ws",0);
  comsolDebug(f_debug," ... done\n");
}

/** This function loads the actual workspace (state of the structure simulation)<br>
    First deletes all variables*/
//void load_checkpoint(clEnv *env,const char *name){
//  comsolDebug(f_debug," Loading checkpoint %s ...\n", name);
//  const char command[50] = "flload comsol_checkpoint.mph";
//  strcpy( &(command[5]) , name);
//  // first we delete everything from the workspace
//  clEvalExpr(env,"clear all",0);
//  // this command should evaluate in the workspace the following line:
//  // "loade checkpointfilename" then a "checkpointfilename.ws" will be loaded in the workspace
//  // the workspace contains all the information about the current state of the simulation
//  clEvalExpr(env,command,0);
//  comsolDebug(f_debug," ... done\n");
//}

// ==================== THE MAIN SIMULATION FUNCTION ===========
CL_EXPORT void work
(
  clEnv* env,
  int nOut,
  clData* out[],
  int nIn,
  clData* in[] )
{
  int i=3;                /* Loop variable */
  int j;                  /* Loop variable */
  int size_coord_entries; /* Size of coordinate entries (x and y) */
  int size_elements;      /* Size of surface mesh elements */

  double* coords;        /* Coordinates of surface mesh nodes */
  double* elems;         /* ??? */
  double* d_x;
  double* d_y;
  double* old_vel_x;     /* Old velocity x components */
  double* old_vel_y;     /* Old velocity y components */
  double* old_displ_x;   /* Old displacement x components */
  double* old_displ_y;   /* Old displacement y components */

  clData* cl_coords;          /* Comsol list of coordinates */
  clData* cl_elems;           /* Comsol list of elements */
  clData* fd_x;
  clData* fd_y;
  clData* cldata_old_vel_x;   /* Comsol list of old velocity x components */
  clData* cldata_old_vel_y;   /* Comsol list of old velocity y components */
  clData* cldata_old_displ_x; /* Comsol list of old displacement x components */
  clData* cldata_old_displ_y; /* Comsol list of old displacement y components */

  clData* in_arg[5];
  clData* out_arg[5];
  clData* it_num;
  FILE* f_x_coords;
  FILE* f_y_coords;
  char filename[50];

  struct hostent* pHostInfo;             /* holds info about a machine */
  struct sockaddr_in Address;            /* Internet socket address stuct */
  long nHostAddress;                       /* long address of the host we want to connect to */
  int nAddressSize = sizeof(struct sockaddr_in);

  int redo_step , doStep , checkpoint_number;
  int succesStep = 1;
  int save_iteration_nr = 0;
  int coresp_measure_point = -1;

  double oldDisplacement[2];
  double displacementDelta[2];
  double oldVelocity[2];
  double velocityDelta[2];

  double* x_coords_with_duplic = NULL;
  double* y_coords_with_duplic = NULL;
  double* x_coords = NULL;
  double* y_coords = NULL;

  int* comsol_to_fsice_node_indices = NULL;
  int* fsice_to_comsol_nodes_indices = NULL;
  int size_nodes_wdupl;
  int size_faces;
  int size_nodes;
  int is_duplicate;
  int size_duplicate_nodes;
  double distance = 0.0;
  double distance_sub1 = 0.0;
  double distance_sub2 = 0.0;

  //delete & create the comsol debug file
  comsolOpen(f_debug, "comsol_debug.tmp", "w");
  //comsolClose(f_debug);
  //comsolOpen(f_debug, "simulation_debug.tmp","w");
  comsolOpen(f_x_coords,"simulation_x_coords.tmp","w");
  comsolOpen(f_y_coords,"simulation_y_coords.tmp","w");

  comsolDebug(f_debug, "Setting up socket connection to preCICE driver...\n");
  socketToPrec = socket(AF_INET,SOCK_STREAM,0);
  if(socketToPrec == SOCKET_ERROR){
    comsolDebug(f_debug,"Could not make a socket\n");
  }
  Address.sin_addr.s_addr=INADDR_ANY;
  Address.sin_port=htons(COMMUNICATION_PORT);
  Address.sin_family=AF_INET;
  comsolDebug(f_debug,"Binding to Prec port %d\n",COMMUNICATION_PORT);
  if(bind( socketToPrec , (struct sockaddr*)&Address,sizeof(Address)) == SOCKET_ERROR) {
    comsolDebug(f_debug,"Could not connect to host\n");
  }
  getsockname( socketToPrec, (struct sockaddr *) &Address,(socklen_t *)&nAddressSize);
  comsolDebug(f_debug,"opened socket as fd (%d) on port (%d) for stream i/o\n",socketToPrec, ntohs(Address.sin_port) );
  if(listen( socketToPrec , QUEUE_SIZE ) == SOCKET_ERROR){
    comsolDebug(f_debug,"Could not listen\n");
  }
  comsolDebug(f_debug,"Waiting for Precice side \n");
  CommSocket = accept( socketToPrec , (struct sockaddr*)&Address , (socklen_t *)&nAddressSize );
  comsolDebug(f_debug,"Connection recieved from Precice \n");
  // =================== END OF SOCKET SETUP ==================

  //clEvalExpr(env, SIM_LOAD, 0); // loading the model , and solve for "0" forces
  // --- load checkpoint ---
  FSI_Recv_int_through_socket(CommSocket , &checkpoint_number);
  if (checkpoint_number > 0){
    comsolDebug(f_debug," Load checkpoint ...\n");
    //sprintf(filename,"comsol_checkpoint");
    //load_checkpoint(env,&(filename[0]));
    clEvalExpr(env,"clear all",0);
    if (clGetLastError(env) != NULL){
      comsolDebug(f_debug, " Error clearing all: %s\n", clGetLastError(env));
    }
    clEvalExpr(env, "flload comsol_checkpoint.mph", 0);
    //clEvalExpr(env, "load comsol_checkpoint.ws", 0);
    if (clGetLastError(env) != NULL){
      comsolDebug(f_debug, " Error loading checkpoint: %s\n", clGetLastError(env));
    }
    clEvalExpr(env, SIM_RESTART, 0);
    if (clGetLastError(env) != NULL){
      comsolDebug(f_debug, " Error initializing restart: %s\n", clGetLastError(env));
    }
    comsolDebug(f_debug, " ... done\n");
  }
  else {
    comsolDebug(f_debug," Creating scenario ...\n");
    clEvalExpr(env, SIM_LOAD, 0);
    if (clGetLastError(env) != NULL){
      comsolDebug(f_debug, " Error creating simulation: %s\n", clGetLastError(env));
    }
    comsolDebug(f_debug, " ... done\n");
  }

  // ============   WE NEED A 1D GRID   ===========
  clEvalExpr(env, "displacement_x = posteval(fem,'u','Refine',2,'Edim',1);", 0);
  if (clGetLastError(env) != NULL){
    comsolDebug(f_debug, " Error evaluating displacement_x: %s\n", clGetLastError(env));
  }

  cl_coords = clEvalExpr(env, "displacement_x.p;", 1);
  if (clGetLastError(env) != NULL){
    comsolDebug(f_debug, " Error fetching displacement x values: %s\n", clGetLastError(env));
  }

  cl_elems = clEvalExpr(env, "displacement_x.t;", 1);

  size_coord_entries = clGetNElems(cl_coords);
  size_elements = clGetNElems(cl_elems);
  // we divide the number of coordinates by 2 and the nr of triangles by 3 (because these lengths show the linear size)
  comsolDebug(f_debug, " Num nodes %d FACES: %d ,Type %d ==== \n",
              (size_coord_entries/2)+1 , size_elements/2 , clGetType(cl_elems) );

  coords = clGetRealPtr(cl_coords);
  elems = clGetRealPtr(cl_elems);

  x_coords_with_duplic = g_malloc((size_coord_entries) * sizeof(double));
  y_coords_with_duplic = g_malloc((size_coord_entries) * sizeof(double));
  x_coords = g_malloc((size_coord_entries) * sizeof(double));
  y_coords = g_malloc((size_coord_entries) * sizeof(double));
  comsol_to_fsice_node_indices = g_malloc((size_coord_entries) * sizeof(int));
  fsice_to_comsol_nodes_indices = g_malloc((size_coord_entries) * sizeof(int));

  size_nodes_wdupl = size_coord_entries / 2;
  size_faces = size_elements / 2;
  comsolDebug(f_debug, "size_coord_entries = %d, size_elements = %d, size_faces = %d\n",
              size_coord_entries, size_elements, size_faces);

  comsolDebug(f_debug, "Assigning node coords...\n");
  for (i=0; i < size_nodes_wdupl; i++){
    comsolDebug(f_debug, "   x = %f, y = %f\n", coords[i*2], coords[(i*2)+1]);
    x_coords_with_duplic[i] = coords[i*2];
    y_coords_with_duplic[i] = coords[(i*2)+1];
  }

  comsolDebug(f_debug, "Setting up node maps Comsol <-> FSI*ce ...\n");

  /* Since the FSIce mesh consists of trianges, it has to use a dummy node, which
     is set at the origin. */
  x_coords[0] = 0.0;
  y_coords[0] = 0.0;
  size_nodes = 1;
  fsice_to_comsol_nodes_indices[0] = -1;
  size_duplicate_nodes = 0;
  for (i=0; i < size_nodes_wdupl; i++){
    is_duplicate = 0;
    comsolDebug(f_debug, "i = %d, unique node count = %d, coords = %f, %f\n",
                i, size_nodes, x_coords_with_duplic[i], y_coords_with_duplic[i]);
    for (j=0; j < size_nodes; j++){ /* Check, if duplicate */
//      distance_sub1 = pow ( x_coords_with_duplic[i] - x_coords[j], 2 );
//      distance_sub2 = pow ( y_coords_with_duplic[i] - y_coords[j], 2 );
//      distance = sqrt ( distance_sub1 + distance_sub2 );
      distance = dist(x_coords_with_duplic[i], x_coords[j],
                      y_coords_with_duplic[i], y_coords[j]);
      if (distance < epsilon){
        comsolDebug (
            f_debug, "Same node as those at index %d with coords %f, %f (= %f, %f)\n",
            j, x_coords[j], y_coords[j],
            x_coords_with_duplic[i], y_coords_with_duplic[i] );
	      comsol_to_fsice_node_indices[i] = j;
	      is_duplicate = 1;
	      size_duplicate_nodes ++;
	      break;
	    }
    }
    if (is_duplicate == 0){ /* If not a duplicate */
      x_coords[size_nodes] = x_coords_with_duplic[i];
      y_coords[size_nodes] = y_coords_with_duplic[i];
      fsice_to_comsol_nodes_indices[size_nodes] = i;
      comsol_to_fsice_node_indices[i] = size_nodes;
      //comsol_to_fsice_node_indices[i] = comsol_to_fsice_node_indices[i] - size_duplicate_nodes;
      size_nodes ++;
    }
  }
  comsolDebug(f_debug, "Comsol -> FSIce node indices:\n");
  for (i=0; i < size_nodes_wdupl; i++){
       comsolDebug(f_debug, "  %d -> %d\n", i, comsol_to_fsice_node_indices[i]);
  }
  comsolDebug(f_debug, "FSIce -> Comsol node indices:\n");
  for (i=0; i < size_nodes; i++){
     comsolDebug(f_debug, "  %d -> %d\n", i, fsice_to_comsol_nodes_indices[i]);
  }

  /* =========== Create mesh and data ============== */
  comsolDebug(f_debug, "Creating mesh nodes...\n");
  mesh = FSI_Mesh_new(size_nodes, size_faces, 0, "COMSOL_MESH");
  /*FSI_Mesh_set_node ( mesh , 0 , 0.0 , 0.0 , 0.0 );*/
  for (i=0; i < size_nodes; i++){ /* Create nodes */
     FSI_Mesh_set_node(mesh, i, x_coords[i], y_coords[i], 0.0);
     comsolDebug(f_debug, "  Mesh node %d coords: %f, %f\n", i, x_coords[i], y_coords[i]);
  }

  comsolDebug(f_debug, "Subtracting 1 from Comsol element node indices...\n");
  for (i=0; i < size_elements; i++){
     comsolDebug(f_debug, "   Subtract 1 from elems[%d] = %d\n", i, (int)elems[i]);
     elems[i]--;
  }

  comsolDebug ( f_debug, "Creating mesh triangles...\n" );
  double x1, x2, y1, y2;
  int index_node_1;
  int index_node_2;
  for ( i=0; i < size_faces; i++ ) { /* Create triangles (also in 2D) */
    comsolDebug ( f_debug, "   elems[%d] = %d, elems[%d] = %d\n",
                  i*2, (int)elems[i*2], (i*2)+1, (int)elems[(i*2)+1] );
    index_node_1 = comsol_to_fsice_node_indices[(int)elems[i*2]];
    index_node_2 = comsol_to_fsice_node_indices[(int)elems[(i*2)+1]];
    comsolDebug ( f_debug, "   Triangle node 1 = %d, node2 = %d\n",
                  index_node_1, index_node_2 );
    x1 = x_coords[index_node_1];
    y1 = y_coords[index_node_1];
    x2 = x_coords[index_node_2];
    y2 = y_coords[index_node_2];
    if ( revertDirection(x1, y1, x2, y2, 0.0, 0.0) > 0 ) {
      FSI_Mesh_triangle_new ( mesh, i, index_node_2, index_node_1, 0 );
      comsolDebug ( f_debug, "   Revert Mesh triangle: %d , %d , %d \n",
                    index_node_2, index_node_1, 0 );
      comsolDebug ( f_debug, "   Original Mesh triangle: %d , %d , %d \n",
                    (int)elems[(i*2)+1], (int)elems[i*2], 0 );
    }
    else {
      FSI_Mesh_triangle_new ( mesh, i, index_node_1, index_node_2, 0 );
      comsolDebug ( f_debug, "   Mesh triangle: %d , %d , %d \n",
                    index_node_1, index_node_2, 0 );
      comsolDebug ( f_debug, "   Original Mesh triangle: %d , %d , %d \n",
                    (int)elems[i*2], (int)elems[(i*2)+1], 0 );
    }
  }

  //FSI_Mesh_save_vtk_FSI_Mesh( mesh ,"Benchmark_noforces.vtk");
  //comsolDebug(f_debug," Mesh created and saved ... \n");

  mesh_forces = FSI_Mesh_Data_new(mesh, 1, 3, "Forces");
  mesh_displacements = FSI_Mesh_Data_new(mesh, 1, 3, "Displacements");
  mesh_displacementdeltas = FSI_Mesh_Data_new(mesh, 1, 3, "DisplacementDeltas");
  mesh_velocities = FSI_Mesh_Data_new(mesh, 1, 3, "Velocities");
  mesh_velocitydeltas = FSI_Mesh_Data_new(mesh, 1, 3, "VelocityDeltas");

  // ======== HERE WE SET EVERYTHIG(Quantities) TO ZERO =======
  for (i=0; i < size_nodes; i++){
    FSI_Data_set_vector(mesh_forces, i, 0.0, 0.0, 0.0);
    FSI_Data_set_vector(mesh_displacements, i, 0.0, 0.0, 0.0);
    FSI_Data_set_vector(mesh_displacementdeltas, i, 0.0, 0.0, 0.0);
    FSI_Data_set_vector(mesh_velocities, i, 0.0, 0.0, 0.0);
    FSI_Data_set_vector(mesh_velocitydeltas, i, 0.0, 0.0, 0.0);
  }

  // send the complete MESH
  comsolDebug(f_debug, " Sending mesh to ComsolPrecice ...\n");
  FSI_Send_mesh_socket(mesh, CommSocket);
  FSI_Send_quantity_socket(mesh, "Displacements", CommSocket);
  FSI_Send_quantity_socket(mesh, "DisplacementDeltas", CommSocket);
  FSI_Send_quantity_socket(mesh, "Velocities", CommSocket);
  FSI_Send_quantity_socket(mesh, "VelocityDeltas", CommSocket);
  comsolDebug(f_debug, " ... done\n");

  // recieve forcesnextstep_simulation();
  //comsolDebug(f_debug," Receiving mesh from ComsolPrecice... \n");
  //FSI_Recv_quantity_socket(mesh, CommSocket);

//  // --- load checkpoint ---
//  FSI_Recv_int_through_socket(CommSocket , &checkpoint_number);
//  if (checkpoint_number > 0){
//    sprintf(filename,"comsol_checkpoint");
//    load_checkpoint(env,&(filename[0]));
//    clEvalExpr(env, SIM_RESTART, 0);
//  }

  // it is important in these stage to call this function !!!
  //comsolDebug(f_debug," Calling nextstep_simulation ... \n");
  nextstep_simulation();
  //clEvalFunc(env,"nextstep_simulation", 0 , out_arg , 0 , in_arg);

# ifdef Debug
  sprintf(filename,"comsol_init.vtk");
  FSI_Mesh_save_vtk_FSI_Mesh(mesh , filename);
# endif

  clEvalExpr(env, SIM_FORCES, 0);

  //comsolDebug(f_debug," Calling nextstep_simulation ... \n");
  nextstep_simulation();
  //clEvalFunc(env,"nextstep_simulation", 0 , out_arg , 0 , in_arg); // reset the pointer counters

/*
 - recieve flag to continue
 yes
 - recieve flag to restart or to make a new step
 - recieve the length of the time step
 - recieve forces
 - simulate (new or repeat step, with the forces and with the given time step)
 - send flag with the status
 - send displacements, velocity
 no
 - shut down
*/
  //======================== MAIN TIMESTEPPING LOOP ===================
  //
  // Overview:
  //
  // receive do
  // while(do)
  //    recveive redo
  //    receive dt
  //    receive fsice mesh
  //    if (not redo)
  //       fetch comsol values (as old)
  //    if (redo)
  //       comsol redo
  //    else
  //       comsol compute next
  //    fetch comsol values (as new)
  //    write comsol values to fsi mesh
  //    receive do checkpoint
  //    if (do checkpoint)
  //       write checkpoint
  //    receive do
  FSI_Recv_int_through_socket(CommSocket , &doStep);
  iteration_number = 0;
  comsolDebug(f_debug, "\n ====== Main timestepping loop ======\n");
  while (doStep > 0){
    comsolDebug(f_debug," Recieving status , it number : %d \n",iteration_number);
    FSI_Recv_int_through_socket(CommSocket , &redo_step);
    FSI_Recv_double_through_socket(CommSocket , &supervisor_timestep);
    comsolDebug(f_debug, " redo_step=%d, supervisor_timestep=%f\n", redo_step, supervisor_timestep);

    comsolDebug(f_debug," Recieving forces ... \n");
    FSI_Recv_quantity_socket(mesh, CommSocket);

    if (!(redo_step > 0)){
      if (iteration_number > 0){
        comsolDebug(f_debug," Free memory from last iteration ...\n");
        g_free(cldata_old_vel_x);
        g_free(cldata_old_vel_y);
        g_free(old_vel_x);
        g_free(old_vel_y);

        g_free(cldata_old_displ_x);
        g_free(cldata_old_displ_y);
        g_free(old_displ_x);
        g_free(old_displ_y);
        comsolDebug(f_debug," ... done\n");
      }
      comsolDebug(f_debug," Retrieve old velocities and displ from comsol ...\n");
      clEvalExpr(env,"velocity_x = posteval(fem,'ut','Refine',2,'Edim',1);",0);
      //comsolDebug(f_debug," 1 ...\n");
      cldata_old_vel_x = clEvalExpr(env,"velocity_x.d;",1);
      //comsolDebug(f_debug," 2 ...\n");
      old_vel_x = clGetRealPtr(cldata_old_vel_x);
      //comsolDebug(f_debug," 3 ...\n");
      clEvalExpr(env,"velocity_y = posteval(fem,'vt','Refine',2,'Edim',1);",0);
      //comsolDebug(f_debug," 4 ...\n");
      cldata_old_vel_y = clEvalExpr(env,"velocity_y.d;",1);
      //comsolDebug(f_debug," 5 ...\n");
      old_vel_y = clGetRealPtr(cldata_old_vel_y);

      clEvalExpr(env,"displacement_x = posteval(fem,'u','Refine',2,'Edim',1);",0);
      cldata_old_displ_x = clEvalExpr(env,"displacement_x.d;",1);
      old_displ_x = clGetRealPtr ( cldata_old_displ_x );
      clEvalExpr(env,"displacement_y = posteval(fem,'v','Refine',2,'Edim',1);",0);
      cldata_old_displ_y = clEvalExpr(env,"displacement_y.d;",1);
      old_displ_y = clGetRealPtr ( cldata_old_displ_y );
      comsolDebug(f_debug," ... done\n");
    }

    /*if (iteration_number >= 50 ) {*/
    if (redo_step > 0){ // repeat the step
      comsolDebug(f_debug," Repeat Step\n");
      clEvalExpr(env , SIM_NEXT_STEP_RESTART , 0);
      comsolDebug(f_debug," ... done\n");
    }
    else { // step forward
      comsolDebug(f_debug," Step forward ...\n");
      clEvalExpr(env , SIM_NEXT_STEP , 0);
      comsolDebug(f_debug," ... done\n ");
    }
    /*}*/
    // to reset the point counters
    comsolDebug(f_debug, " Compute timestep ...\n");
    nextstep_simulation();
    //clEvalFunc(env,"nextstep_simulation", 0 , out_arg , 0 , in_arg);
    comsolDebug(f_debug," ... done\n");
    comsolDebug(f_debug," Get computed displacements ...\n");
    // ============ Evaluating the displacements and the diff of displacements ====
    // get the displacement in X direction
    clEvalExpr(env,"displacement_x = posteval(fem,'u','Refine',2,'Edim',1);",0);
    fd_x = clEvalExpr(env,"displacement_x.d;",1);
    d_x = clGetRealPtr(fd_x);
    // get the displacement in Y direction
    clEvalExpr(env,"displacement_y = posteval(fem,'v','Refine',2,'Edim',1);",0);
    fd_y = clEvalExpr(env,"displacement_y.d;",1);
    d_y = clGetRealPtr(fd_y);
    comsolDebug(f_debug," ... done\n");
//    for ( i=0 ; i<size_nodes_wdupl ; i++ ){ // if we take always the difference then it is good for implicit/explicit
//      xDisplold = FSI_Data_get_value( mesh_displacement, comsol_to_fsice_node_indices[i+1] , 1 );
//      yDisplold = FSI_Data_get_value( mesh_displacement, comsol_to_fsice_node_indices[i+1] , 2 );
//      //comsolDebug(f_debug,"   Old displ = %f, %f \n", xDisplold, yDisplold );
//      xDisplold = d_x[i] - xDisplold;
//      yDisplold = d_y[i] - yDisplold;
//      //comsolDebug(f_debug,"   New displ = %f, %f \n", xDisplold, yDisplold );
//      FSI_Data_set_vector( mesh_displacement , comsol_to_fsice_node_indices[i+1] , d_x[i] , d_y[i] , 0 );
//	   FSI_Data_set_vector( mesh_diffdisplacement , comsol_to_fsice_node_indices[i+1] , xDisplold , yDisplold , 0 );
//      //double x = FSI_Mesh_get_node_x (mesh, comsol_to_fsice_node_indices[i+1]);
//      //double y = FSI_Mesh_get_node_y (mesh, comsol_to_fsice_node_indices[i+1]);
//      //comsolDebug(f_debug,"   Coordinates = %f, %f \n\n", x, y);
//    }
    comsolDebug(f_debug, " Set displ in FSIce mesh ...\n");
    for (i=1; i < size_nodes; i++){ // if we take always the difference then it is good for implicit/explicit
      oldDisplacement[0] = old_displ_x[fsice_to_comsol_nodes_indices[i]];
      oldDisplacement[1] = old_displ_y[fsice_to_comsol_nodes_indices[i]];
//      xDisplold = FSI_Data_get_value ( mesh_displacements, i+1 , 1 );
//      yDisplold = FSI_Data_get_value ( mesh_displacements, i+1 , 2 );
      //comsolDebug(f_debug,"   Old displ = %f, %f \n", xDisplold, yDisplold );
      displacementDelta[0] = d_x[fsice_to_comsol_nodes_indices[i]] - oldDisplacement[0];
      displacementDelta[1] = d_y[fsice_to_comsol_nodes_indices[i]] - oldDisplacement[1];
      comsolDebug(f_debug,"Set displacement %d = %f, %f \n", i,
                  d_x[fsice_to_comsol_nodes_indices[i]],
                  d_y[fsice_to_comsol_nodes_indices[i]]);
      FSI_Data_set_vector(mesh_displacements, i,
                          d_x[fsice_to_comsol_nodes_indices[i]],
                          d_y[fsice_to_comsol_nodes_indices[i]] , 0);
      FSI_Data_set_vector(mesh_displacementdeltas, i,
                          displacementDelta[0], displacementDelta[1], 0);
      //double x = FSI_Mesh_get_node_x (mesh, comsol_to_fsice_node_indices[i+1]);
      //double y = FSI_Mesh_get_node_y (mesh, comsol_to_fsice_node_indices[i+1]);
      //comsolDebug(f_debug,"   Coordinates = %f, %f \n\n", x, y);
    }
    comsolDebug(f_debug," ... done\n");

     // free the unsused memory ?!
     g_free(fd_x);
     g_free(d_x);
     g_free(fd_y);
     g_free(d_y);

     // ============ Saving the velocoties ======================
     comsolDebug(f_debug," Get computed velocities ... \n");
     clEvalExpr(env,"velocity_x = posteval(fem,'ut','Refine',2,'Edim',1);",0);
     fd_x = clEvalExpr(env,"velocity_x.d;",1);
     d_x = clGetRealPtr(fd_x);
     // get the velocities in Y direction
     clEvalExpr(env,"velocity_y = posteval(fem,'vt','Refine',2,'Edim',1);",0);
     fd_y = clEvalExpr(env,"velocity_y.d;",1);
     d_y = clGetRealPtr(fd_y);
     comsolDebug(f_debug," ... done \n");
     comsolDebug(f_debug, " Set velocities in FSIce mesh ...\n");
     for ( i=1; i < size_nodes; i++ ) {
       oldVelocity[0] = old_vel_x[fsice_to_comsol_nodes_indices[i]];
       oldVelocity[1] = old_vel_y[fsice_to_comsol_nodes_indices[i]];
       //comsolDebug(f_debug,"   Old velocity = %f, %f \n", oldVelocity[0], oldVelocity[1] );
       //comsolDebug(f_debug,"   New velocity = %f, %f \n", d_x[fsice_to_comsol_nodes_indices[i]], d_y[fsice_to_comsol_nodes_indices[i]] );
       oldVelocity[2] = 0.0;
       velocityDelta[0] = d_x[fsice_to_comsol_nodes_indices[i]] - oldVelocity[0];
       velocityDelta[1] = d_y[fsice_to_comsol_nodes_indices[i]] - oldVelocity[1];
       //comsolDebug(f_debug,"   Velocity delta = %f, %f \n", velocityDelta[0], velocityDelta[1] );
       FSI_Data_set_vector(mesh_velocities, i,
                           d_x[fsice_to_comsol_nodes_indices[i]],
                           d_y[fsice_to_comsol_nodes_indices[i]],
                           0.0);
       FSI_Data_set_vector(mesh_velocitydeltas, i, velocityDelta[0],
                           velocityDelta[1], 0.0);
//       FSI_Data_set_vector( mesh_velocitydeltas, i+1, 0.0, 0.0, 0.0 );
       double x = FSI_Mesh_get_node_x(mesh, i);
       double y = FSI_Mesh_get_node_y(mesh, i);
       //comsolDebug(f_debug,"   Coordinates = %f, %f \n\n", x, y);
     }

//     for (i=0 ; i<size_nodes_wdupl; i++ ) {
//       oldVelocity[0] = FSI_Data_get_value (mesh_velocity, comsol_to_fsice_node_indices[i+1], 1);
//       oldVelocity[1] = FSI_Data_get_value (mesh_velocity, comsol_to_fsice_node_indices[i+1], 2);
//       comsolDebug(f_debug,"   Old velocity = %f, %f \n", oldVelocity[0], oldVelocity[1] );
//       comsolDebug(f_debug,"   New velocity = %f, %f \n", d_x[i], d_y[i] );
//       oldVelocity[2] = 0.0;
//       velocityDelta[0] = d_x[i] - oldVelocity[0];
//       velocityDelta[1] = d_y[i] - oldVelocity[1];
//       comsolDebug(f_debug,"   Velocity delta = %f, %f \n", velocityDelta[0], velocityDelta[1] );
//       FSI_Data_set_vector( mesh_velocity, comsol_to_fsice_node_indices[i+1], d_x[i], d_y[i], 0.0 );
//       //FSI_Data_set_vector( mesh_velocitydelta, comsol_to_fsice_node_indices[i+1], velocityDelta[0],
//       //                     velocityDelta[1], 0.0 );
//       FSI_Data_set_vector( mesh_velocitydelta, comsol_to_fsice_node_indices[i+1], 0.0,
//                            0.0, 0.0 );
//       double x = FSI_Mesh_get_node_x (mesh, comsol_to_fsice_node_indices[i+1]);
//       double y = FSI_Mesh_get_node_y (mesh, comsol_to_fsice_node_indices[i+1]);
//       comsolDebug(f_debug,"   Coordinates = %f, %f \n\n", x, y);
//     }

     comsolDebug(f_debug," Send Status ... \n");
     FSI_Send_int_through_socket(CommSocket , &succesStep);
     comsolDebug(f_debug," Send Quantities ... \n");
     FSI_Send_quantity_socket(mesh, "Displacements", CommSocket);
     FSI_Send_quantity_socket(mesh, "DisplacementDeltas", CommSocket);
     FSI_Send_quantity_socket(mesh, "Velocities", CommSocket);
     FSI_Send_quantity_socket(mesh, "VelocityDeltas", CommSocket);

     // Save the actual grid
//#    ifdef Debug
//     sprintf(filename,"Benchmark_COM_%d.vtk" , iteration_number);
//     FSI_Mesh_save_vtk_FSI_Mesh( mesh , filename );
//#    endif

     // ========  save checkpoint  ===========
     FSI_Recv_int_through_socket(CommSocket , &checkpoint_number);
     if (checkpoint_number > 0){
     	 sprintf(filename,"comsol_checkpoint");
     	 save_checkpoint(env,&(filename[0]));
     }

     // Recieve/Check if we move to the next simulation step
     FSI_Recv_int_through_socket(CommSocket , &doStep);
     iteration_number++;

     // free the unsused memory ?!
     g_free(fd_x);
     g_free(d_x);
     g_free(fd_y);
     g_free(d_y);
  }
  // ================================ SIMULATION ENDS HERE =============

  close(CommSocket);
  comsolDebug(f_debug," After simulation cycle ... doing cleanup\n");

  // geomobject*/
  //to clean up the variable on the comsol side
  clEvalFunc (env, "cleanup_simulation", 0, out_arg, 0, in_arg);

  if (mesh != NULL) FSI_Mesh_destroy(mesh);

  comsolDebug(f_debug , " End of simulation \n");
  comsolClose(f_debug);
  fclose(f_x_coords);
  fclose(f_y_coords);
}
