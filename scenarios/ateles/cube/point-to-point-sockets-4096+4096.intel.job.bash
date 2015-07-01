#!/bin/bash

#@ job_name = ateles/cube
#@ job_type = MPICH

#@ initialdir = $(home)

###@ output = job.out
###@ error  = job.err

#@ wall_clock_limit = 00:20:00

#@ network.MPI = sn_all,not_shared,us

#@ class          = general
#@ island_count   = 1,2
#@ node           = 512
#@ tasks_per_node = 16
#@ queue

 distribution_type='point-to-point'
communication_type='sockets'

export configuration="${distribution_type}-${communication_type}"
export       LEFT_NP='4096'
export      RIGHT_NP='4096'

SCENARIO_DIR="${HOME}/scenarios/ateles/cube"
     RUN_DIR="${HOME}/runs/ateles/cube"
   # RUN_DIR="${WORK}/runs/ateles/cube"

. "${SCENARIO_DIR}/ateles.intel.bash"
