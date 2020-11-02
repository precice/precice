PROGRAM main
  IMPLICIT NONE
  CHARACTER*512                   :: config
  CHARACTER*50                    :: participantName, meshName, writeInitialData, readItCheckp, writeItCheckp
  CHARACTER*50                    :: readDataName, writeDataName
  INTEGER                         :: rank, commsize, ongoing, dimensions, meshID, bool, numberOfVertices, i,j
  INTEGER                         :: readDataID, writeDataID
  DOUBLE PRECISION                :: dt
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: vertices, writeData, readData
  INTEGER, DIMENSION(:), ALLOCATABLE :: vertexIDs

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
    writeDataName = 'dataOne'
    readDataName = 'dataTwo'
  ENDIF
  IF(participantName .eq. 'SolverTwo') THEN
    writeDataName = 'dataTwo'
    readDataName = 'dataOne'
  ENDIF

  rank = 0
  commsize = 1
  dt = 1
  numberOfVertices = 3
  CALL precicef_create(participantName, config, rank, commsize)

  ! Allocate dummy mesh with only one vertex
  CALL precicef_get_dims(dimensions)
  ALLOCATE(vertices(numberOfVertices*dimensions))
  ALLOCATE(vertexIDs(numberOfVertices))
  ALLOCATE(readData(numberOfVertices*dimensions))
  ALLOCATE(writeData(numberOfVertices*dimensions))
  CALL precicef_get_mesh_id(meshName, meshID)

  do i = 1,numberOfVertices,1
    do j = 1,dimensions,1
      vertices((i - 1)*dimensions + j ) = i-1
      readData((i - 1)*dimensions + j ) = i-1
      writeData((i - 1)*dimensions + j ) = i-1
    enddo
    vertexIDs(i) = i-1
  enddo

  CALL precicef_set_vertices(meshID, numberOfVertices, vertices, vertexIDs)
  DEALLOCATE(vertices)

  CALL precicef_get_data_id(readDataName,meshID,readDataID)
  CALL precicef_get_data_id(writeDataName,meshID,writeDataID)

  CALL precicef_initialize(dt)

  CALL precicef_is_action_required(writeInitialData, bool)
  IF (bool.EQ.1) THEN
    WRITE (*,*) 'DUMMY: Writing initial data'
  ENDIF
  CALL precicef_initialize_data()

  CALL precicef_is_coupling_ongoing(ongoing)
  DO WHILE (ongoing.NE.0)
  
    CALL precicef_is_action_required(writeItCheckp, bool)
    
    IF (bool.EQ.1) THEN
      WRITE (*,*) 'DUMMY: Writing iteration checkpoint'
      CALL precicef_mark_action_fulfilled(writeItCheckp)
    ENDIF

    CALL precicef_is_read_data_available(bool)
    IF (bool.EQ.1) THEN
      CALL precicef_read_bvdata(readDataID, numberOfVertices, vertexIDs, readData)
    ENDIF

    WRITE (*,*) 'readData: ', readData

    writeData = readData + 1

    CALL precicef_is_write_data_required(dt, bool)
    IF (bool.EQ.1) THEN
      CALL precicef_write_bvdata(writeDataID, numberOfVertices, vertexIDs, writeData)
    ENDIF

    CALL precicef_advance(dt)

    CALL precicef_is_action_required(readItCheckp, bool)
    IF (bool.EQ.1) THEN
      WRITE (*,*) 'DUMMY: Reading iteration checkpoint'
      CALL precicef_mark_action_fulfilled(readItCheckp)
    ELSE
      WRITE (*,*) 'DUMMY: Advancing in time'
    ENDIF

    CALL precicef_is_coupling_ongoing(ongoing)

  ENDDO

  CALL precicef_finalize()
  WRITE (*,*) 'DUMMY: Closing Fortran solver dummy...'

  DEALLOCATE(writeData)
  DEALLOCATE(readData)
  DEALLOCATE(vertexIDs)

END PROGRAM
