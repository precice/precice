export FSI_FSIDUMMYA_JAVA=off
export FSI_FSIDUMMYA_WORKERS=20
export FSI_FSIDUMMYA_XML=./FSIDummies/workspace/FSIDummyWorkbench.xml
export FSI_FSIDUMMYA_BUFFER_SIZE=4096
mpiexec -np 2 ./dummyA  2 1

