#include "precice/adapters/c/SolverInterfaceC.h"
#include "precice/adapters/c/Constants.h"
#include <stdio.h>
#include <stdlib.h>
#include <rpc/rpc.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <unistd.h>

#include "fsi_mesh.h"
#include "fsi_mesh_utils.h"
#include "fsi_map-mat.h"
#include "fsi_interface_socket.h"

// on this socket is the communication between COMSOL and Precise trough FSI
#define COMMUNICATION_SOCKET       53219
#define HOSTNAME                  "localhost"
#define MAX_NR_POINTS              10000

static int socketToComsol;
static FSI_Mesh * mesh = NULL;
static FSI_Data * meshForces = NULL;
static FSI_Data * meshDisplacements = NULL;
static FSI_Data * meshDisplacementDeltas = NULL;
static FSI_Data * meshVelocities = NULL;
static FSI_Data * meshVelocityDeltas = NULL;

static int iterationNumber;
static int maxIterationNumber = 20;
static double timeStep = 0.01;
static int redoStep = 0;
static int doStep = 1;
static int succesStep = 1;
static int forcesID = -1;
static int velocitiesID = -1;
static int displacementsID = -1;
static int displacementDeltasID = -1;
static int velocityDeltasID = -1;
static int meshID = -1;
static int dimensions = 0;

//static double x_actualshift[MAX_NR_POINTS];
//static double y_actualshift[MAX_NR_POINTS];
static FILE *fd;


void readForces
(
  FSI_Mesh * mesh,
  FSI_Data * data,
  int        dataID )
{
  int i;
  double coords[dimensions];
  double value[dimensions];

  fprintf ( fd, "Reading forces (dim = %d)\n", dimensions );

  if ( dimensions == 2 ) {
    for ( i=1; i < FSI_Mesh_get_num_nodes(mesh); i++ ) {
      coords[0] = FSI_Mesh_get_node_x ( mesh, i );
      coords[1] = FSI_Mesh_get_node_y ( mesh, i );
      precicec_readVectorData ( dataID, i-1, value );
      FSI_Data_set_vector ( data, i, value[0], value[1], 0.0 );
    }
  }
  else if ( dimensions == 3 ) {
    for ( i=1; i < FSI_Mesh_get_num_nodes(mesh); i++ ) {
      coords[0] = FSI_Mesh_get_node_x ( mesh , i );
      coords[1] = FSI_Mesh_get_node_y ( mesh, i );
      coords[2] = FSI_Mesh_get_node_z ( mesh, i );
      precicec_readVectorData ( dataID, i-1, value );
      FSI_Data_set_vector ( data, i, value[0], value[1], value[2] );
    }
  }
}

void writePreciceData
(
  FSI_Mesh * mesh,
  FSI_Data * data,
  int        dataID )
{
  int i;
  int iDim;
  double coords[dimensions];
  double value[dimensions];

  fprintf ( fd, "Writing data \"%s\"\n", data->_name );

  if ( dimensions == 2 ) {
    for ( i=1; i < FSI_Mesh_get_num_nodes(mesh); i++ ) {
      coords[0] = FSI_Mesh_get_node_x ( mesh, i );
      coords[1] = FSI_Mesh_get_node_y ( mesh, i );
      value[0] = FSI_Data_get_value ( data, i, 1 );
      value[1] = FSI_Data_get_value ( data, i, 2 );
      fprintf ( fd, "Set value %d = %e,%e at coords = %e, %e\n",
                i, value[0], value[1], coords[0], coords[1] );
      precicec_writeVectorData ( dataID, i-1, value );
    }
  }
  else if ( dimensions == 3 ) {
    for ( i=1; i < FSI_Mesh_get_num_nodes(mesh); i++ ) {
      coords[0] = FSI_Mesh_get_node_x ( mesh, i );
      coords[1] = FSI_Mesh_get_node_y ( mesh, i );
      coords[2] = FSI_Mesh_get_node_z ( mesh, i );
      value[0] = FSI_Data_get_value ( data, i, 1 );
      value[1] = FSI_Data_get_value ( data, i, 2 );
      value[2] = FSI_Data_get_value ( data, i, 3 );
      fprintf ( fd, "Set value %d = %e,%e, %e at coords = %e, %e, %e\n",
                i, value[0], value[1], value[2], coords[0], coords[1], coords[2] );
      precicec_writeVectorData ( dataID, i-1, value );
    }
  }
}

//void writeVelocities
//(
//  FSI_Mesh * mesh,
//  FSI_Data * data )
//{
//  int i;
//  double coords[2]
//  double value[2];
//
//  for ( i=1; i < FSI_Mesh_get_num_nodes(mesh_tmp); i++ ) {
//    coords[0] = FSI_Mesh_get_node_x ( mesh, i );
//    coords[1] = FSI_Mesh_get_node_y ( mesh, i );
//    value[0] = FSI_Data_get_value ( data, i, 1 );
//    value[1] = FSI_Data_get_value ( data, i, 2 );
//    fprintf ( fd, "Set value(velocity): %d : %e,%e with: %e,%e \n",
//              i, coords[0], coords[1], value[0], value[1] );
//    precicec_writeVectorData ( velocityID, coords, value );
//  }
//}
//
//void writeVelocityDeltas
//(
//  FSI_Mesh * mesh,
//  FSI_Data * data )
//{
//   int i;
//   double coords[2];
//   double velDelta[2];
//
//   for ( i=1; i < FSI_Mesh_get_num_nodes(mesh); i++ ) {
//     coords[0] = FSI_Mesh_get_node_x ( mesh, i );
//     coords[1] = FSI_Mesh_get_node_y ( mesh, i );
//     velDelta[0] = FSI_Data_get_value ( data, i, 1 );
//     velDelta[1] = FSI_Data_get_value ( data, i, 2 );
//     fprintf ( fd,"Setting velocity delta: %d : %e,%e with: %e,%e \n",
//               i, coords[0], coords[1], velDelta[0], velDelta[1] );
//     precicec_writeVectorData ( velocityDeltaID, coords, velDelta );
//   }
//}
//
//void writeDisplacements(FSI_Mesh * mesh_tmp , FSI_Data * mesh_data_new_diff , FSI_Data * mesh_data_displacement ){
//  int i;
//  double coords[2] , displacement[2];
//  for (i = 1 ; i < FSI_Mesh_get_num_nodes( mesh_tmp ) ; i++) {
//    coords[0] = FSI_Mesh_get_node_x( mesh_tmp , i );
//    coords[1] = FSI_Mesh_get_node_y( mesh_tmp , i );
////    displacement[0] = FSI_Data_get_value( mesh_data_displacement , i , 1 );
////    displacement[1] = FSI_Data_get_value( mesh_data_displacement , i , 2 );
//    displacement[0] = FSI_Data_get_value( mesh_data_new_diff , i , 1 );
//    displacement[1] = FSI_Data_get_value( mesh_data_new_diff , i , 2 );
//    fprintf(fd,"Setting displacements: %d : %e,%e with: %e,%e \n",
//            i ,coords[0] , coords[1] , displacement[0] , displacement[1] );
//    precicec_writeVectorData ( displacementID, coords , displacement );
//  }
//}
//
//void writeDisplacementDeltas(FSI_Mesh * mesh_tmp , FSI_Data * mesh_data_new_diff , FSI_Data * mesh_data_displacement ){
//  int i;
//  double coords[2] , displacement[2];
//  for (i = 1 ; i < FSI_Mesh_get_num_nodes( mesh_tmp ) ; i++) {
//    coords[0] = FSI_Mesh_get_node_x( mesh_tmp , i );
//    coords[1] = FSI_Mesh_get_node_y( mesh_tmp , i );
////    displacement[0] = FSI_Data_get_value( mesh_data_displacement , i , 1 );
////    displacement[1] = FSI_Data_get_value( mesh_data_displacement , i , 2 );
//    displacement[0] = FSI_Data_get_value( mesh_data_new_diff , i , 1 );
//    displacement[1] = FSI_Data_get_value( mesh_data_new_diff , i , 2 );
//    fprintf(fd,"Setting displacements: %d : %e,%e with: %e,%e \n",
//            i ,coords[0] , coords[1] , displacement[0] , displacement[1] );
//    precicec_writeVectorData ( displacementID, coords , displacement );
//  }
//}

//====================================== MAIN =========================================
int main(int argc, char** argv)
{
  int hServerSocket, hServerSocket1; /* handle to socket */
  struct hostent * pHostInfo;        /* holds info about a machine */
  struct sockaddr_in Address;       /* Internet socket address stuct */
  long nHostAddress;       /* long address of the host we want to connect to */
  int fHostPort, sHostPort;
  int nAddressSize = sizeof(struct sockaddr_in);
	double coords[2], value[2];
//	int pointNr , trinagleNr;
	int i;
  int checkpointNumber;
	char filename[50];

	fd = fopen("comsolprecice.tmp", "w");

  if (argc < 1){
    // fprintf(fd,"TOO FEW PARAMETERS: \n -accessorName [IN] Name of the solver accessing the interface. Has to match one of the names specified in the configuration xml file\n");
    fprintf(fd,"TOO FEW PARAMETERS: \n");
    fprintf(fd," -configFileName [IN] (Path and) name of the xml configuration file containing the precice configuration.");
	  return 0;
	}

  fprintf(fd, "Calling precicec_createSolverInterface with: Comsol, %s \n", argv[1]);
  precicec_createSolverInterface("Comsol", argv[1], 0, 1);
  dimensions = precicec_getDimensions();

  socketToComsol = socket(AF_INET, SOCK_STREAM, 0);
  pHostInfo = gethostbyname(HOSTNAME); /* get IP address from name */
  /* copy address into long */
  memcpy(&nHostAddress, pHostInfo->h_addr, pHostInfo->h_length);
  /* fill address struct */
  Address.sin_addr.s_addr = nHostAddress;
  Address.sin_port = htons(COMMUNICATION_SOCKET);
  Address.sin_family = AF_INET; /* AF_INET represents the address family INET for Internet sockets. */
  fprintf(fd,"Connecting to %d on port %d \n", nHostAddress, sHostPort);
  /* connect to host */
  printf("Connecting to Comsol...\n"); fflush(stdout);
  while(connect(socketToComsol, (struct sockaddr*) &Address, sizeof(Address)) == -1){
    usleep(1000000);
  }
	printf("Connection with success ... Now recieving Mesh\n"); fflush(stdout);

  // --- send the checkpoint number which should be loaded ---
  checkpointNumber = -1; // -1 means it no checkpoints will be loaded
                          // only when checkpoint_number is greater than 0,
                          // will a checkpoint loaded
  if (precicec_isActionRequired(precicec_actionReadSimulationCheckpoint())){
    checkpointNumber = 1;
    precicec_fulfilledAction(precicec_actionReadSimulationCheckpoint());
  }
  FSI_Send_int_through_socket(socketToComsol, &checkpointNumber);

  /* GET THE MESH FROM COMSOL AND MAKE INITIALIZATIONS */
	mesh = FSI_Mesh_new_empty(); /* this initialization is needed !!!! */
	fprintf(fd, "Receiving FSIce mesh ... \n");
	FSI_Recv_mesh_socket(mesh , socketToComsol);

  fprintf(fd, "Building Precice mesh ... \n");
  int meshID = precicec_getMeshID("WetSurface");
//  precicec_createMeshBuilder ( "Geometry" );
	/* create the vertices in the Precice geometry */
	for (i=1; i < FSI_Mesh_get_num_nodes(mesh); i++){
    coords[0] = FSI_Mesh_get_node_x(mesh, i);
    coords[1] = FSI_Mesh_get_node_y(mesh, i);
    fprintf(fd, "Adding vertex with coodrs = %f, %f\n", coords[0], coords[1]);
    int index = precicec_setMeshVertex(meshID, coords);
    fprintf(fd, "   ... returned index = %d\n", index);
	}
	/* Create the edges in the geometry */
	for (i=0; i < FSI_Mesh_get_num_triangles(mesh); i++){
	  fprintf(fd, "Adding edge with vertex indices = %d, %d\n",
	          FSI_Mesh_get_face_node1(mesh,i)-1, FSI_Mesh_get_face_node2(mesh,i)-1);
	  precicec_setMeshEdge(meshID,
                        FSI_Mesh_get_face_node1(mesh,i)-1,
                        FSI_Mesh_get_face_node2(mesh,i)-1);
	}

  fprintf(fd, "Recieve and create data ...\n");
  meshForces = FSI_Mesh_Data_new(mesh, 1, 3, "Forces");
	FSI_Recv_quantity_socket(mesh, socketToComsol);
  FSI_Recv_quantity_socket(mesh, socketToComsol);
  FSI_Recv_quantity_socket(mesh, socketToComsol);
  FSI_Recv_quantity_socket(mesh, socketToComsol);
	meshDisplacements = FSI_Mesh_get_Data(mesh, "Displacements");
  meshDisplacementDeltas = FSI_Mesh_get_Data(mesh, "DisplacementDeltas");
  meshVelocities = FSI_Mesh_get_Data(mesh, "Velocities");
  meshVelocityDeltas = FSI_Mesh_get_Data(mesh, "VelocityDeltas");

  //fprintf(fd, "Sending mesh to comsol ...\n");
  //FSI_Send_quantity_socket(mesh, "Forces", socketToComsol);

  printf("Initializing Coupling ... \n"); fflush(stdout);
	timeStep = precicec_initialize();

	forcesID = precicec_getDataID("Forces");
	if (precicec_hasData("Velocities")){
	  velocitiesID = precicec_getDataID("Velocities");
	}
	if (precicec_hasData("VelocityDeltas")){
	  velocityDeltasID = precicec_getDataID("VelocityDeltas");
	}
	if (precicec_hasData("Displacements")){
	  displacementsID = precicec_getDataID("Displacements");
	}
	if (precicec_hasData("DisplacementDeltas")){
	  displacementDeltasID = precicec_getDataID("DisplacementDeltas");
	}
  fprintf(fd, "Exporting precice geometry ...\n");
  precicec_exportMesh("Precice");

  // MAIN TIMESTEPPING LOOP
  //
  // Overview:
  //
  // if (do)
  //    send do to comsol
  //    send redo to comsol
  //    send dt to comsol
  //    copy forces from precice to fsice mesh
  //    send forces to comsol
  //    receive comsol values
  //    copy comsol values from fsice to precice mesh
  //    send do checkpoint to comsol
  // send do to comsol
  // send do checkpoint to comsol
  //
  while (precicec_isCouplingOngoing() > 0){
    iterationNumber++;
	  if (precicec_isActionRequired("write-iteration-checkpoint")){
      precicec_fulfilledAction("write-iteration-checkpoint");
	  }
	  if (precicec_isActionRequired("read-iteration-checkpoint")){
		  redoStep = 1; /* repeat the timestep */
	    precicec_fulfilledAction("read-iteration-checkpoint");
	  }
	  else {
	    redoStep = 0; /* do not repeat the timestep */
	  }
	  fprintf(fd, "Sending coupling state to Comsol ...\n");
    FSI_Send_int_through_socket(socketToComsol, &doStep);
    FSI_Send_int_through_socket(socketToComsol, &redoStep);
    FSI_Send_double_through_socket(socketToComsol, &timeStep);

    readForces(mesh, meshForces, forcesID);
    fprintf(fd, "Sending force data to Comsol ...\n");
    FSI_Send_quantity_socket(mesh, "Forces", socketToComsol);

	  fprintf(fd, "Receiving Comsol data 1 ...\n");
	  FSI_Recv_int_through_socket(socketToComsol , &succesStep);
	  fprintf(fd, "Receiving Comsol data 2 ...\n");
	  FSI_Recv_quantity_socket(mesh, socketToComsol);
    fprintf(fd, "Receiving Comsol data 3 ...\n");
    FSI_Recv_quantity_socket(mesh, socketToComsol);
    fprintf(fd, "Receiving Comsol data 4 ...\n");
    FSI_Recv_quantity_socket(mesh, socketToComsol);
    fprintf(fd, "Receiving Comsol data 5 ...\n");
    FSI_Recv_quantity_socket(mesh, socketToComsol);
    fprintf(fd, "... done\n");

    if (velocitiesID != -1){
      writePreciceData(mesh, meshVelocities, velocitiesID);
    }
    if (velocityDeltasID != -1){
      writePreciceData(mesh, meshVelocityDeltas, velocityDeltasID);
    }
    if (displacementsID != -1){
	    writePreciceData(mesh, meshDisplacements, displacementsID);
    }
    if (displacementDeltasID != -1){
	    writePreciceData(mesh, meshDisplacementDeltas, displacementDeltasID);
    }

	  // ------- Exchange data -----------
    fprintf ( fd, "Calling preCICE advance ... \n");
	  timeStep = precicec_advance(timeStep);

    // --- send the checkpoint number where this should be stored ---
    // if the checkpoint < 1 then no checkpoint will be written
    checkpointNumber = -1;
    if (precicec_isActionRequired(precicec_actionWriteSimulationCheckpoint())){
      checkpointNumber = 1;
      precicec_fulfilledAction(precicec_actionWriteSimulationCheckpoint());
    }
    FSI_Send_int_through_socket(socketToComsol , &checkpointNumber);
  }
	fprintf(fd,"Tell COMSOL simulation ended  ... \n");
  // --- send the checkpoint number where this should be stored ---
  // if the checkpoint < 1 then no checkpoint will be written
//  checkpointNumber = -1;
//  if (precicec_isActionRequired(precicec_actionWriteSimulationCheckpoint())){
//    checkpointNumber = 1;
//    precicec_fulfilledAction(precicec_actionWriteSimulationCheckpoint());
//  }
//  FSI_Send_int_through_socket(socketToComsol , &checkpointNumber);

	doStep = 0;
  FSI_Send_int_through_socket ( socketToComsol , &doStep );
  close ( socketToComsol );
	fprintf ( fd, "End simulation ... \n" );
	fclose ( fd );
  return 1;
}
