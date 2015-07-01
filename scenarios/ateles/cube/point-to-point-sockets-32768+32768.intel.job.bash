#!/bin/bash

#@ job_name = ateles/cube
#@ job_type = MPICH

#@ initialdir = $(home)/bin/cube

###@ output = job.out
###@ error  = job.err

#@ wall_clock_limit = 00:30:00

#@ network.MPI = sn_all,not_shared,us

#@ class          = special
#@ island_count   = 8,16
#@ node           = 4096
#@ tasks_per_node = 16
#@ queue

 distribution_type='point-to-point'
communication_type='sockets'

export configuration="${distribution_type}-${communication_type}"
export       LEFT_NP='32768'
export      RIGHT_NP='32768'

SCENARIO_DIR="${HOME}/scenarios/ateles/cube"
     RUN_DIR="${HOME}/runs/ateles/cube"
   # RUN_DIR="${WORK}/runs/ateles/cube"

. "${SCENARIO_DIR}/ateles.intel.bash"
