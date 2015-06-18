#!/bin/bash

export ATELES="${SCENARIO_DIR}/ateles.ibm"

mkdir -p                                           "${RUN_DIR}"
cp    -f -r "${SCENARIO_DIR}/left"                 "${RUN_DIR}"
cp    -f -r "${SCENARIO_DIR}/right"                "${RUN_DIR}"
cp    -f    "${SCENARIO_DIR}/${configuration}.xml" "${RUN_DIR}/precice-config.xml"

cd "${RUN_DIR}"

rm -r -f .*.address*

chmod +x  'left/ateles.ibm.subjob.bash'
chmod +x 'right/ateles.ibm.subjob.bash'

cat << EOF > 'cmdfile'
 left/ateles.ibm.subjob.bash@1%${LEFT_NP}%mpich2:*
right/ateles.ibm.subjob.bash@2%${RIGHT_NP}%mpich2:*
wait
complete
EOF

. '/etc/profile'
. '/etc/profile.d/modules.sh'

module unload mpi.intel
module   load mpi.ibm
module   load lrztools

export MP_NEWJOB=parallel
export MP_LABELIO=yes
export MP_PGMMODEL=mpmd

poe -cmdfile 'cmdfile'
