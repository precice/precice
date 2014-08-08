export FSI_FSIDUMMYA_JAVA=off
export FSI_FSIDUMMYA_WORKERS=20
export FSI_FSIDUMMYA_XML=/work/atanasoa/gitlab/precice/tools/communication_dummies/FSIDummies/workspace/FSIDummyWorkbench.xml
export FSI_FSIDUMMYA_BUFFER_SIZE=4096
/home/atanasoa/intel/impi/4.1.3.049/intel64/bin/mpiexec -np 6 ./dummyA  2 3

