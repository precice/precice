      PROGRAM main
      IMPLICIT NONE
      CHARACTER*100 name
      CHARACTER*100 config
      CHARACTER*(*) meshName
      PARAMETER (meshName = 'WetSurface')
      CHARACTER*50 writeInitialData
      CHARACTER*50 readSimCheckp
      CHARACTER*50 writeSimCheckp
      CHARACTER*50 readItCheckp
      CHARACTER*50 writeItCheckp
      INTEGER rank
      INTEGER commsize
      INTEGER ongoing
      INTEGER meshID
      INTEGER bool
      DOUBLE PRECISION dtlimit
            
      CALL precicef_action_write_initial_data(writeInitialData)
      CALL precicef_action_read_sim_checkp(readSimCheckp)
      CALL precicef_action_write_sim_checkp(writeSimCheckp)
      CALL precicef_action_read_iter_checkp(readItCheckp)
      CALL precicef_action_write_iter_checkp(writeItCheckp)
      
      WRITE (*,*) 'Starting Fortran solver dummy...'
      CALL getarg(1,name)
      CALL getarg(2,config)
      rank = 0
      commsize = 1
      CALL precicef_create(name, config, rank, commsize)      
      CALL precicef_initialize(dtlimit)            

      CALL precicef_action_required(writeInitialData, bool)
      IF (bool.EQ.1) THEN
        WRITE (*,*) 'Writing initial data'
        CALL precicef_initialize_data()
      ENDIF

      CALL precicef_action_required(readSimCheckp, bool)
      IF (bool.EQ.1) THEN
        WRITE (*,*) 'Reading simulation checkpoint'
        CALL precicef_fulfilled_action(readSimCheckp)
      ENDIF

      CALL precicef_ongoing(ongoing)
      DO WHILE (ongoing.NE.0)
        WRITE (*,*) 'Action name:', writeItCheckp, '|- end'
        CALL precicef_action_required(writeItCheckp, bool)
        IF (bool.EQ.1) THEN
          WRITE (*,*) 'Writing iteration checkpoint'
          CALL precicef_fulfilled_action(writeItCheckp)
        ENDIF
      
		    CALL precicef_advance(dtlimit)
        CALL precicef_ongoing(ongoing)
        
        CALL precicef_action_required(writeSimCheckp, bool)
        IF (bool.EQ.1) THEN
          WRITE (*,*) 'Writing simulation checkpoint'
          CALL precicef_fulfilled_action(writeSimCheckp)
        ENDIF
        
        CALL precicef_action_required(readItCheckp, bool)
        IF (bool.EQ.1) THEN
          WRITE (*,*) 'Reading iteration checkpoint'
          CALL precicef_fulfilled_action(readItCheckp)
        ELSE
          WRITE (*,*) 'Advancing in time'
        ENDIF
      ENDDO
      CALL precicef_finalize()
      
      END PROGRAM 
