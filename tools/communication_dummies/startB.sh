export FSI_FSIDUMMYB_JAVA=off
export FSI_FSIDUMMYB_WORKERS=20
export FSI_FSIDUMMYB_XML=/work/atanasoa/gitlab/precice/tools/communication_dummies/FSIDummies/workspace/FSIDummyWorkbench.xml
export FSI_FSIDUMMYB_BUFFER_SIZE=4096
/home/atanasoa/intel/impi/4.1.3.049/intel64/bin/mpiexec -np 6 ./dummyB  2 3
