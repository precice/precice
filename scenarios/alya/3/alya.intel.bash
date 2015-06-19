#!/bin/bash

export ALYA="${SCENARIO_DIR}/alya.x.intel"

mkdir -p                                           "${RUN_DIR}"
cp    -f -r "${SCENARIO_DIR}/nastin"               "${RUN_DIR}"
cp    -f -r "${SCENARIO_DIR}/solidz"               "${RUN_DIR}"
cp    -f    "${SCENARIO_DIR}/${configuration}.xml" "${RUN_DIR}/precice-config.xml"

cd "${RUN_DIR}"

rm -r -f .*.address*

L1='1'
L2="${NASTIN_NP}"

sed -n -e "${L1},${L2}p" "${LOADL_HOSTFILE}" > 'nastin/hostfile'

L1="$((NASTIN_NP + 1        ))"
L2="$((NASTIN_NP + SOLIDZ_NP))"

sed -n -e "${L1},${L2}p" "${LOADL_HOSTFILE}" > 'solidz/hostfile'

. '/etc/profile'
. '/etc/profile.d/modules.sh'

module unload mpi.ibm
module   load mpi.intel/4.1

export SUBJOB

pushd 'nastin'
  mpiexec -n "${NASTIN_NP}" -f 'hostfile' "${ALYA}" 'c' &> "${configuration}-$((NASTIN_NP + SOLIDZ_NP))-${NASTIN_NP}.log" &
popd

pushd 'solidz'
  mpiexec -n "${SOLIDZ_NP}" -f 'hostfile' "${ALYA}" 'c' &> "${configuration}-$((NASTIN_NP + SOLIDZ_NP))-${SOLIDZ_NP}.log"
popd
