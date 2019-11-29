PROGRAM main
  IMPLICIT NONE
  CHARACTER*50                    :: config, solverParameters, participantName, meshName, writeInitialData 
  CHARACTER*50                    :: DataOneName, DataTwoName
  CHARACTER*50                    :: readItCheckp, writeItCheckp
  INTEGER                         :: rank, commsize, ongoing, dimensions, meshID, bool, n, i,j
  INTEGER			  :: DataOne_ID, DataTwo_ID
  REAL                            :: dtlimit
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: vertex, DataOne, DataTwo
  INTEGER, DIMENSION(:), ALLOCATABLE :: vertexID
      
  ! Constants in f90 have to be prefilled with blanks to be compatible with preCICE
  writeInitialData(1:50)='                                                  '
  readItCheckp(1:50)='                                                  '
  writeItCheckp(1:50)='                                                  '
  
  CALL precicef_action_write_initial_data(writeInitialData)
  CALL precicef_action_read_iter_checkp(readItCheckp)
  CALL precicef_action_write_iter_checkp(writeItCheckp)

  WRITE (*,*) 'DUMMY: Starting Fortran solver dummy...'
  CALL getarg(1, config)
  CALL getarg(2, participantName)
  CALL getarg(3, meshName)

  IF(participantName .eq. 'SolverOne') THEN  
    DataOneName = 'Forces'
    DataTwoName = 'Velocities'
  ENDIF
  IF(participantName .eq. 'SolverTwo') THEN  
    DataOneName = 'Velocities'
    DataTwoName = 'Forces'
  ENDIF

  rank = 0
  commsize = 1
  dtlimit = 1
  n = 3
  CALL precicef_create(participantName, config, rank, commsize)

  ! Allocate dummy mesh with only three vertices 
  CALL precicef_get_dims(dimensions)
  ALLOCATE(vertex(n*dimensions))
  ALLOCATE(vertexID(n))
  vertex = 0
  vertexID = 0

  CALL precicef_get_mesh_id(meshName, meshID)
  print*, 'meshID = ', meshID

  do i = 1,n,1
    do j = 1,dimensions,1
      vertex((i - 1)*dimensions + j ) = i
    enddo
    vertexID(i) = i
  enddo

  CALL precicef_set_vertices(meshID, n*dimensions, vertex, vertexID)  
  DEALLOCATE(vertex)     

  CALL precicef_get_data_id(DataOneName,meshID,DataOne_ID)
  CALL precicef_get_data_id(DataTwoName,meshID,DataTwo_ID)

  ALLOCATE(DataOne(n))
  ALLOCATE(DataTwo(n))
  DataOne = 0
  DataTwo = 0

  CALL precicef_initialize(dtlimit)            

  CALL precicef_action_required(writeInitialData, bool)
  IF (bool.EQ.1) THEN
  WRITE (*,*) 'DUMMY: Writing initial data'
  ENDIF
  CALL precicef_initialize_data()

  CALL precicef_ongoing(ongoing)
  DO WHILE (ongoing.NE.0)

    WRITE (*,*) 'DUMMY: Reading iteration checkpoint'
    do i = 1,n,1
      CALL precicef_read_sdata(DataTwo_ID, vertexID(i), DataTwo(i))
    enddo

    DataOne = DataTwo + 1

    WRITE (*,*) 'DUMMY: Writing iteration checkpoint'
    do i=1,n,1
      CALL precicef_write_sdata(DataOne_ID, vertexID(i), DataOne(i));
    enddo

    WRITE (*,*) 'DUMMY: Advancing in time'
    CALL precicef_advance(dtlimit)
    CALL precicef_ongoing(ongoing)

  ENDDO

  CALL precicef_finalize()
  WRITE (*,*) 'DUMMY: Closing Fortran solver dummy...'

END PROGRAM 
