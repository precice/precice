PROGRAM main
  IMPLICIT NONE
  CHARACTER*512                   :: config
  CHARACTER*50                    :: participantName, meshName
  CHARACTER*50                    :: readDataName, writeDataName
  INTEGER                         :: rank, commsize, ongoing, dimensions, bool, numberOfVertices, i,j
  DOUBLE PRECISION                :: dt
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: vertices, writeData, readData
  INTEGER, DIMENSION(:), ALLOCATABLE :: vertexIDs

  WRITE (*,*) 'DUMMY: Starting Fortran solver dummy...'

  CALL getarg(1, config)
  CALL getarg(2, participantName)

  IF(participantName .eq. 'SolverOne') THEN
    writeDataName = 'Data-One'
    readDataName = 'Data-Two'
    meshName = 'SolverOne-Mesh'
  ENDIF
  IF(participantName .eq. 'SolverTwo') THEN
    writeDataName = 'Data-Two'
    readDataName = 'Data-One'
    meshName = 'SolverTwo-Mesh'
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

  do i = 1,numberOfVertices,1
    do j = 1,dimensions,1
      vertices((i - 1)*dimensions + j ) = i-1
      readData((i - 1)*dimensions + j ) = i-1
      writeData((i - 1)*dimensions + j ) = i-1
    enddo
    vertexIDs(i) = i-1
  enddo

  CALL precicef_set_vertices(meshName, numberOfVertices, vertices, vertexIDs)
  DEALLOCATE(vertices)

  CALL precicef_requires_initial_data(bool)
  IF (bool.EQ.1) THEN
    WRITE (*,*) 'DUMMY: Writing initial data'
  ENDIF
  CALL precicef_initialize(dt)

  CALL precicef_is_coupling_ongoing(ongoing)
  DO WHILE (ongoing.NE.0)

    CALL precicef_requires_writing_checkpoint(bool)

    IF (bool.EQ.1) THEN
      WRITE (*,*) 'DUMMY: Writing iteration checkpoint'
    ENDIF

    CALL precicef_read_bvdata(meshName, readDataName, numberOfVertices, vertexIDs, readData)

    WRITE (*,*) 'readData: ', readData

    writeData = readData + 1

    CALL precicef_write_bvdata(meshName, writeDataName, numberOfVertices, vertexIDs, writeData)

    CALL precicef_advance(dt)

    CALL precicef_requires_reading_checkpoint(bool)
    IF (bool.EQ.1) THEN
      WRITE (*,*) 'DUMMY: Reading iteration checkpoint'
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
