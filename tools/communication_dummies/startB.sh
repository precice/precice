export FSI_FSIDUMMYB_JAVA=off
export FSI_FSIDUMMYB_WORKERS=20
export FSI_FSIDUMMYB_XML=./FSIDummies/workspace/FSIDummyWorkbench.xml
export FSI_FSIDUMMYB_BUFFER_SIZE=4096
/home/uekerman/Software/openmpi-1.8.1/bin/mpiexec -np 6 ./dummyB  2 3
