/* Copyright (C) 2011 Technische Universitaet Muenchen
 * This file is part of the preCICE project. For conditions of distribution and
 * use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License */
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
