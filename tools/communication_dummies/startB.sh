export FSI_FSIDUMMYB_JAVA=off
export FSI_FSIDUMMYB_WORKERS=20
export FSI_FSIDUMMYB_XML=./FSIDummies/workspace/FSIDummyWorkbench.xml
export FSI_FSIDUMMYB_BUFFER_SIZE=4096
mpiexec -np 2 ./dummyB  2 1
