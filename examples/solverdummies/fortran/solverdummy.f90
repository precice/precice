PROGRAM main
  IMPLICIT NONE
  CHARACTER*512                   :: config
  CHARACTER*50                    :: participantName, meshName, writeInitialData, readItCheckp, writeItCheckp
  CHARACTER*50                    :: DataReadName, DataWriteName
  INTEGER                         :: rank, commsize, ongoing, dimensions, meshID, bool, N, i,j
  INTEGER                         :: DataRead_ID, DataWrite_ID
  REAL                            :: dtlimit
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: vertex, DataWrite, DataRead
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
    DataWriteName = 'Forces'
    DataReadName = 'Velocities'
  ENDIF
  IF(participantName .eq. 'SolverTwo') THEN  
    DataWriteName = 'Velocities'
    DataReadName = 'Forces'
  ENDIF

  rank = 0
  commsize = 1
  dtlimit = 1
  N = 3             !Number of vertices
  CALL precicef_create(participantName, config, rank, commsize)

  ! Allocate dummy mesh with only one vertex 
  CALL precicef_get_dims(dimensions)
  ALLOCATE(vertex(N*dimensions))
  ALLOCATE(vertexID(N))
  ALLOCATE(DataRead(N*dimensions))
  ALLOCATE(DataWrite(N*dimensions))
  CALL precicef_get_mesh_id(meshName, meshID)

  do i = 1,N,1
    do j = 1,dimensions,1
      vertex((i - 1)*dimensions + j ) = i-1
      DataRead((i - 1)*dimensions + j ) = i-1
      DataWrite((i - 1)*dimensions + j ) = i-1
    enddo
    vertexID(i) = i-1
  enddo

  CALL precicef_set_vertices(meshID, N, vertex, vertexID)  
  DEALLOCATE(vertex)
  
  CALL precicef_get_data_id(DataReadName,meshID,DataRead_ID)
  CALL precicef_get_data_id(DataWriteName,meshID,DataWrite_ID)
        
  CALL precicef_initialize(dtlimit)            

  CALL precicef_action_required(writeInitialData, bool)
  IF (bool.EQ.1) THEN
    WRITE (*,*) 'DUMMY: Writing initial data'
  ENDIF
  CALL precicef_initialize_data()

  CALL precicef_ongoing(ongoing)
  DO WHILE (ongoing.NE.0)
  
    CALL precicef_action_required(writeItCheckp, bool)
    IF (bool.EQ.1) THEN
      CALL precicef_write_bvdata(DataWrite_ID, N, vertexID, DataWrite)
      WRITE (*,*) 'DUMMY: Writing iteration checkpoint'
      CALL precicef_mark_action_fulfilled(writeItCheckp)
    ENDIF

    CALL precicef_advance(dtlimit)
    CALL precicef_ongoing(ongoing)

    CALL precicef_action_required(readItCheckp, bool)
    IF (bool.EQ.1) THEN
      CALL precicef_read_bvdata(DataRead_ID, N, vertexID, DataRead)
      WRITE (*,*) 'DUMMY: Reading iteration checkpoint'
      CALL precicef_mark_action_fulfilled(readItCheckp)
    ELSE
      WRITE (*,*) 'DUMMY: Advancing in time'
    ENDIF

    WRITE (*,*) 'DataRead: ', DataRead

    DataWrite = DataRead + 1
    
  ENDDO
  
  CALL precicef_finalize()
  WRITE (*,*) 'DUMMY: Closing Fortran solver dummy...'

END PROGRAM 
