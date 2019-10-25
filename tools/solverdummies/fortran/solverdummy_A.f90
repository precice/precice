PROGRAM main
  IMPLICIT NONE
  CHARACTER*50                    :: config, participantName, meshName, writeInitialData, readItCheckp, writeItCheckp
  INTEGER                         :: rank, commsize, ongoing, dimensions, meshID, bool, n, i,j
  INTEGER			  :: solver_One_Data_ID, solver_Two_Data_ID
  REAL                            :: dtlimit
  REAL, DIMENSION(:), ALLOCATABLE :: vertex, solver_One_Data, solver_Two_Data
  INTEGER, DIMENSION(:), ALLOCATABLE :: vertexID
      
  ! Constants in f90 have to be prefilled with blanks to be compatible with preCICE
  writeInitialData(1:50)='                                                  '
  readItCheckp(1:50)='                                                  '
  writeItCheckp(1:50)='                                                  '
  
  CALL precicef_action_write_initial_data(writeInitialData)
  CALL precicef_action_read_iter_checkp(readItCheckp)
  CALL precicef_action_write_iter_checkp(writeItCheckp)

  WRITE (*,*) 'DUMMY: Starting Fortran solver dummy...'

  rank = 0
  commsize = 1
  CALL precicef_create('dummy_Solver_One', 'precice-config.xml', rank, commsize)

  WRITE (*,*) 'DUMMY: Starting Fortran solver dummy...'

  CALL precicef_get_dims(dimensions)
  ALLOCATE(vertex(n*dimensions))
  ALLOCATE(vertexID(n))
  vertex = 0
  vertexID = 0
  n = 3
  dtlimit = 1

  CALL precicef_get_mesh_id('dummy_Mesh_One', meshID)

  do i = 1,n,1
    do j = 1,dimensions,1
      vertex((i-1)*dimensions + i ) = i
      print*,'vertex = ', vertex((i-1)*dimensions + i )
    enddo
    vertexID(i) = i - 1
    print*,'vertexID = ', vertexID(i)
  enddo

  CALL precicef_set_vertices(meshID, n, vertex, vertexID)  
  DEALLOCATE(vertex)    

  CALL precicef_get_data_id('solver_One_Data',meshID,solver_One_Data_ID)
  CALL precicef_get_data_id('solver_Two_Data',meshID,solver_Two_Data_ID)
  print*, 'mesh_One ID = ', solver_One_Data_ID, ' and mesh_Two = ', solver_Two_Data_ID

  ALLOCATE(solver_One_Data(n))
  ALLOCATE(solver_Two_Data(n))
  solver_One_Data = 0
  solver_Two_Data = 0

  CALL precicef_initialize(dtlimit)            

  CALL precicef_action_required(writeInitialData, bool)
  IF (bool.EQ.1) THEN
    WRITE (*,*) 'DUMMY: Writing initial data'
  ENDIF
  CALL precicef_initialize_data()

  CALL precicef_ongoing(ongoing)
  DO WHILE (ongoing.NE.0)

    WRITE (*,*) 'DUMMY: Advancing in time'
    do i = 1,n,1
      CALL precicef_read_sdata(solver_Two_Data_ID, vertexID(i), solver_Two_Data(i))
      print*,'solver_Two_Data = ', solver_Two_Data(i)
      print*,'vertexID = ', vertexID(i)
    enddo

    solver_One_Data = solver_Two_Data + 1

    WRITE (*,*) 'DUMMY: Writing iteration checkpoint'
    do i=1,n,1
      print*,'solver_One_Data = ', solver_One_Data(i)
      print*,'vertexID = ', vertexID(i)
      CALL precicef_write_sdata(solver_One_Data_ID, vertexID(i), solver_One_Data(i));
    enddo

    read(*,*)
    CALL precicef_advance(dtlimit)
    CALL precicef_ongoing(ongoing)

  ENDDO
  
  CALL precicef_finalize()
  WRITE (*,*) 'DUMMY: Closing Fortran solver dummy...'

END PROGRAM 
