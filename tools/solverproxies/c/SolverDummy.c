#include "CouplingInterfaceC.h"

#include <stdio.h>

int main ( int argc, char **argv )
{
   double dt = 0.0;

   if (argc != 3) {
      printf ("Usage: ./exec dummy-name config-file-name \n");
      return 1;
   }
   printf ("Running solver dummy with name %s and config file %s\n",
           argv[1], argv[2]);

   precice_createCouplingInterface (argv[1], argv[2]);

   dt = precice_initialize ();

   while ( precice_isCoupledSimulationOngoing() ) {
      dt = precice_advance ( dt );
   }

   precice_finalize ();

   return 1;
}
