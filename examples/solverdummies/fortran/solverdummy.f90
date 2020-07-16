PROGRAM main
  IMPLICIT NONE
  CHARACTER*512                   :: config
  CHARACTER*50                    :: participantName, meshName, writeInitialData, readItCheckp, writeItCheckp
  INTEGER                         :: rank, commsize, ongoing, dimensions, meshID, vertexID, bool
  REAL                            :: dtlimit
  REAL, DIMENSION(:), ALLOCATABLE :: vertex
      
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

  rank = 0
  commsize = 1
  CALL precicef_create(participantName, config, rank, commsize)

  ! Allocate dummy mesh with only one vertex 
  CALL precicef_get_dims(dimensions)
  ALLOCATE(vertex(dimensions))
  vertex = 0
  CALL precicef_get_mesh_id(meshName, meshID)
  CALL precicef_set_vertex(meshID, vertex, vertexID)  
  DEALLOCATE(vertex)    
        
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
      WRITE (*,*) 'DUMMY: Writing iteration checkpoint'
      CALL precicef_mark_action_fulfilled(writeItCheckp)
    ENDIF

    CALL precicef_advance(dtlimit)
    CALL precicef_ongoing(ongoing)

    CALL precicef_action_required(readItCheckp, bool)
    IF (bool.EQ.1) THEN
      WRITE (*,*) 'DUMMY: Reading iteration checkpoint'
      CALL precicef_mark_action_fulfilled(readItCheckp)
    ELSE
      WRITE (*,*) 'DUMMY: Advancing in time'
    ENDIF
    
  ENDDO
  
  CALL precicef_finalize()
  WRITE (*,*) 'DUMMY: Closing Fortran solver dummy...'

END PROGRAM 
