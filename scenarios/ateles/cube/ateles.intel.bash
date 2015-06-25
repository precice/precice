#!/bin/bash

export ATELES="${SCENARIO_DIR}/ateles.intel"

mkdir -p                                           "${RUN_DIR}"
cp    -f -r "${SCENARIO_DIR}/left"                 "${RUN_DIR}"
cp    -f -r "${SCENARIO_DIR}/right"                "${RUN_DIR}"
cp    -f    "${SCENARIO_DIR}/${configuration}.xml" "${RUN_DIR}/precice-config.xml"

cd "${RUN_DIR}"

rm -r -f .*.address*

L1='1'
L2="${LEFT_NP}"

sed -n -e "${L1},${L2}p" "${LOADL_HOSTFILE}" >  'left/hostfile'

L1="$((LEFT_NP + 1       ))"
L2="$((LEFT_NP + RIGHT_NP))"

sed -n -e "${L1},${L2}p" "${LOADL_HOSTFILE}" > 'right/hostfile'

. '/etc/profile'
. '/etc/profile.d/modules.sh'

module unload mpi.ibm
module   load mpi.intel

export SUBJOB

pushd 'left'
  mpiexec -n "${LEFT_NP}"  -f 'hostfile' "${ATELES}" 'ateles_left.lua'  &> "${configuration}-$((LEFT_NP + RIGHT_NP))-${LEFT_NP}.log"  &
popd

pushd 'right'
  mpiexec -n "${RIGHT_NP}" -f 'hostfile' "${ATELES}" 'ateles_right.lua' &> "${configuration}-$((LEFT_NP + RIGHT_NP))-${RIGHT_NP}.log"
popd
