PROGRAM main
  use PreCICE_solver_if_module
  IMPLICIT NONE
  
  ! We need the length of the strings, set this to a meaningful value in your code.
  ! Here assumed that length = 50 (arbitrary).
  CHARACTER*50                    :: config, participantName, meshName, writeInitialData, readItCheckp, writeItCheckp
  INTEGER                         :: rank, commsize, ongoing, dimensions, meshID, vertexID, bool
  REAL(8)                         :: dtlimit
  REAL(8), DIMENSION(:), ALLOCATABLE :: vertex
  integer(kind=c_int)             :: c_accessorNameLength
  integer(kind=c_int)             :: c_configFileNameLength

  ! Constants in f90 have to be prefilled with blanks to be compatible with preCICE
  writeInitialData(1:50)='                                                  '
  readItCheckp(1:50)='                                                  '
  writeItCheckp(1:50)='                                                  '
  
  CALL precicef_action_write_initial_data(writeInitialData, 50)
  CALL precicef_action_read_iter_checkp(readItCheckp, 50)
  CALL precicef_action_write_iter_checkp(writeItCheckp, 50)

  WRITE (*,*) 'DUMMY: Starting Fortran solver dummy...'
  CALL getarg(1, config)
  CALL getarg(2, participantName)
  CALL getarg(3, meshName)

  rank = 0
  commsize = 1
  CALL precicef_create(participantName, config, rank, commsize, 50, 50)

  ! Allocate dummy mesh with only one vertex 
  CALL precicef_get_dims(dimensions)
  ALLOCATE(vertex(dimensions))
  vertex = 0
  CALL precicef_get_mesh_id(meshName, meshID, 50)
  CALL precicef_set_vertex(meshID, vertex, vertexID)  
  DEALLOCATE(vertex)    
        
  CALL precicef_initialize(dtlimit)            

  CALL precicef_action_required(writeInitialData, bool, 50)
  IF (bool.EQ.1) THEN
    WRITE (*,*) 'DUMMY: Writing initial data'
  ENDIF
  CALL precicef_initialize_data()

  CALL precicef_ongoing(ongoing)
  DO WHILE (ongoing.NE.0)
  
    CALL precicef_action_required(writeItCheckp, bool, 50)
    IF (bool.EQ.1) THEN
      WRITE (*,*) 'DUMMY: Writing iteration checkpoint'
      CALL precicef_fulfilled_action(writeItCheckp, 50)
    ENDIF

    CALL precicef_advance(dtlimit)
    CALL precicef_ongoing(ongoing)

    CALL precicef_action_required(readItCheckp, bool, 50)
    IF (bool.EQ.1) THEN
      WRITE (*,*) 'DUMMY: Reading iteration checkpoint'
      CALL precicef_fulfilled_action(readItCheckp, 50)
    ELSE
      WRITE (*,*) 'DUMMY: Advancing in time'
    ENDIF
    
  ENDDO
  
  CALL precicef_finalize()
  WRITE (*,*) 'DUMMY: Closing Fortran solver dummy...'

END PROGRAM 
