#include <iostream>
#include <sstream>
#include "precice/SolverInterface.hpp"
#include "utils/Dimensions.hpp"
#include "utils/GeometryComputations.hpp"

/**
 * @brief For printing to the command line.
 *
 * @param message  An input stream such as: "Blabla = " << variable << ", ..."
 */
#ifdef STRUCTURE_DEBUG_MODE
#define STRUCTURE_DEBUG(message) \
   { \
      std::ostringstream conv; \
      conv << message; \
      std::cout << conv.str() << std::endl; \
   }
#else
#define STRUCTURE_DEBUG(message)
#endif

#define STRUCTURE_INFO(message) \
   { \
      std::ostringstream conv; \
      conv << message; \
      std::cout << conv.str() << std::endl; \
   }

using precice::utils::DynVector;
using namespace precice;
using tarch::la::raw;

/**
 * @brief Runs structure0815
 */
int main ( int argc, char **argv )
{
   STRUCTURE_INFO ( "Starting Structure0815" );

   if ( argc != 2 ) {
      STRUCTURE_INFO ( "Usage: ./structure0815 configurationFileName" );
      abort ();
   }
   std::string configFileName ( argv[1] );
   std::string meshName ( "WetSurface" );

   SolverInterface cplInterface ( "Structure0815", 0, 1 );
   cplInterface.configure ( configFileName );
   double dt = cplInterface.initialize ();

   std::string writeCheckpoint ( constants::actionWriteIterationCheckpoint() );
   std::string readCheckpoint ( constants::actionReadIterationCheckpoint() );

   double computedTime = 0.0;
   int computedTimeSteps = 0;
   int plotNumber = 0;
   int dimensions = cplInterface.getDimensions();
   if (dimensions != 2){
     STRUCTURE_INFO ( "Structure0815 is only implemented in 2D yet!" );
     exit ( -1 );
   }

   double density = 10.0;
   DynVector gravity (dimensions, 0.0); //-9.81;

   STRUCTURE_INFO ( "Density = " << density );
   STRUCTURE_INFO ( "Gravity = " << gravity );

   DynVector translVelocityChange ( dimensions, 0.0 );

   if ( not cplInterface.hasMesh(meshName) ){
     STRUCTURE_INFO ( "Mesh \"" << meshName << "\" required for coupling!" );
     exit ( -1 );
   }
   MeshHandle handle = cplInterface.getMeshHandle ( meshName );

   int velocitiesID = -1;
   int forcesID = -1;
   int displacementsID = -1;
   int displDeltasID = -1;
   int velocityDeltasID = -1;

   if ( ! cplInterface.hasData(constants::dataForces()) ) {
     STRUCTURE_INFO ( "Data \"Forces\" required for coupling!" );
     exit ( -1 );
   }
   forcesID = cplInterface.getDataID ( constants::dataForces() );
   if ( cplInterface.hasData(constants::dataVelocities()) ) {
     velocitiesID = cplInterface.getDataID ( constants::dataVelocities() );
   }
   if ( cplInterface.hasData(constants::dataDisplacements()) ) {
     displacementsID = cplInterface.getDataID ( constants::dataDisplacements() );
   }
   if ( cplInterface.hasData("DisplacementDeltas") ) {
     displDeltasID = cplInterface.getDataID ( "DisplacementDeltas" );
   }
   if ( cplInterface.hasData("VelocityDeltas") ) {
     velocityDeltasID = cplInterface.getDataID ( "VelocityDeltas" );
   }

   DynVector zero(dimensions, 0.0);
   std::vector<DynVector> velocities ( handle.vertices().size(), zero );
   std::vector<DynVector> velocityChanges ( handle.vertices().size(), zero );
   std::vector<DynVector> displacements ( handle.vertices().size(), zero );
   std::vector<DynVector> displacementChanges ( handle.vertices().size(), zero );

   while ( cplInterface.isCouplingOngoing() ){
      if ( cplInterface.isActionRequired(writeCheckpoint) ){
         cplInterface.fulfilledAction ( writeCheckpoint );
      }
      //if ( computedTimeSteps >= 0 ) { // Start computing deformations
      DynVector coords0 ( zero );
      DynVector coords1 ( zero );
      DynVector centerOfGravity ( zero );
      double totalMass = 0.0;
      double area = 0.0;
      double totalArea = 0.0;
      size_t weights = 0;
      EdgeHandle edges = handle.edges();
      foriter ( EdgeIterator, it, edges ){
        for (int i=0; i<dimensions; i++) coords0[i] = it.vertexCoords(0)[i];
        for (int i=0; i<dimensions; i++) coords1[i] = it.vertexCoords(1)[i];
        typedef precice::utils::GeometryComputations GeoComp;
        area = GeoComp::triangleArea(zero, coords0, coords1);
        area = std::abs ( area ); // since it comes out signed from cross-prod
        totalArea += area;
        if ( not tarch::la::equals(area, 0.0) ) {
          weights++;
          centerOfGravity += (coords0 + coords1) * area / 3.0;
        }
      }
      totalMass = totalArea * density;
      centerOfGravity /= totalArea;

      DynVector force ( zero );
      DynVector totalForce ( gravity );
      totalForce *= totalMass; // Makes gravity acceleration a force
      DynVector r ( zero );
      double torque = 0.0;
      VertexHandle vertices = handle.vertices ();
      for ( VertexIterator it = vertices.begin(); it != vertices.end(); it++ ){
        cplInterface.readVectorData ( forcesID, it.vertexID(), raw(force) );
        totalForce += force;
        for (int i=0; i<dimensions; i++)
          r[i] = -1.0 * centerOfGravity[i] + it.vertexCoords()[i];
        //r = -1.0 * centerOfGravity + dwrap(it.vertexCoords());
        torque += (r[0] * force[1]) - (r[1] * force[0]);
      }

      STRUCTURE_DEBUG ( "Total mass = " << totalMass );
      STRUCTURE_DEBUG ( "Total force = " << totalForce );
      STRUCTURE_DEBUG ( "Center of gravity = " << centerOfGravity );
      STRUCTURE_DEBUG ( "Torque = " << torque );

      // Compute values of next timestep
      translVelocityChange = (totalForce / totalMass) * dt;

      // Set values of next timestep
      DynVector rotForce ( zero );
      DynVector rotVelocityChange ( zero );
      DynVector writeVelocity ( zero );
      DynVector writeDisplacement ( zero );
      double normR = 0.0;
      size_t iVertex = 0;
      foriter ( VertexIterator, it, vertices ){
        // Compute and write velocity
        for (int i=0; i<dimensions; i++)
          r[i] = -1.0 * centerOfGravity[i] + it.vertexCoords()[i];
        normR = tarch::la::norm2 ( r );
        rotForce[0] = torque * -1.0 * r[1];
        rotForce[1] = torque * r[0];
        rotVelocityChange = (rotForce / totalMass) * dt;
        velocityChanges[iVertex] = rotVelocityChange + translVelocityChange;
        writeVelocity = velocities[iVertex] + velocityChanges[iVertex];
        displacementChanges[iVertex] =
            (velocities[iVertex] + velocityChanges[iVertex] * dt * 0.5) * dt;
        writeDisplacement = displacements[iVertex] + displacementChanges[iVertex];

        if ( velocitiesID != -1 ) {
          cplInterface.writeVectorData ( velocitiesID, it.vertexID(), raw(writeVelocity) );
        }
        if ( displacementsID != -1 ) {
          cplInterface.writeVectorData ( displacementsID, it.vertexID(), raw(writeDisplacement) );
        }
        if ( displDeltasID != -1 ) {
          cplInterface.writeVectorData ( displDeltasID, it.vertexID(), raw(displacementChanges[iVertex]) );
        }
        if ( velocityDeltasID != -1 ) {
          cplInterface.writeVectorData ( velocityDeltasID, it.vertexID(), raw(velocityChanges[iVertex]) );
        }
        iVertex++;
      }
      std::ostringstream stream;
      stream << "structure0815-regular." << plotNumber;
      cplInterface.exportMesh ( stream.str() );
      plotNumber ++;

      dt = cplInterface.advance ( dt );

      if ( cplInterface.isActionRequired(readCheckpoint) ) {
         STRUCTURE_DEBUG ( "Loading checkpoint" );
         cplInterface.fulfilledAction ( readCheckpoint );
      }
      else {
         STRUCTURE_DEBUG ( "Adding value changes to values." );
         for ( size_t i=0; i < velocities.size(); i++ ) {
            velocities[i] += velocityChanges[i];
            displacements[i] += displacementChanges[i];
         }
         computedTime += dt;
         computedTimeSteps ++;
      }

      STRUCTURE_DEBUG ( "Computed time = " << computedTime
                        << ", computed timesteps = " << computedTimeSteps );
   }

   STRUCTURE_DEBUG ( "Finalizing coupling" );
   cplInterface.finalize ();

   STRUCTURE_INFO ( "Exiting Structure0815" );

   return 1;
}
