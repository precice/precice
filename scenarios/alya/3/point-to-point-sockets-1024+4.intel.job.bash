#!/bin/bash

#@ job_name = alya/3
#@ job_type = MPICH

#@ initialdir = $(home)

###@ output = job.out
###@ error  = job.err

#@ wall_clock_limit = 01:00:00

#@ network.MPI = sn_all,not_shared,us

#@ class          = general
#@ island_count   = 1
#@ node           = 65
#@ tasks_per_node = 16
#@ queue

 distribution_type='point-to-point'
communication_type='sockets'

export configuration="${distribution_type}-${communication_type}"
export     NASTIN_NP='1024'
export     SOLIDZ_NP='4'

SCENARIO_DIR="${HOME}/scenarios/alya/3"
     RUN_DIR="${HOME}/runs/alya/3"
   # RUN_DIR="${WORK}/runs/alya/3"

. "${SCENARIO_DIR}/alya.intel.bash"
