module PreCICE_solver_if_module
  use, intrinsic :: iso_c_binding
  implicit none

  interface

    subroutine precicef_create(participantName, configFileName, &
      &                        solverProcessIndex, solverProcessSize, &
      &                        lengthAccessorName, lengthConfigFileName) &
      &  bind(c, name='precicef_create_')

      use, intrinsic :: iso_c_binding
      character(kind=c_char), dimension(*) :: participantName
      character(kind=c_char), dimension(*) :: configFileName
      integer(kind=c_int) :: solverProcessIndex
      integer(kind=c_int) :: solverProcessSize
      integer(kind=c_int), value :: lengthAccessorName
      integer(kind=c_int), value :: lengthConfigFileName
    end subroutine precicef_create


    subroutine precicef_initialize(timestepLengthLimit) &
      &  bind(c, name='precicef_initialize_')

      use, intrinsic :: iso_c_binding
      real(kind=c_double) :: timestepLengthLimit
    end subroutine precicef_initialize

    subroutine precicef_advance(timestepLengthLimit) &
      &  bind(c, name='precicef_advance_')

      use, intrinsic :: iso_c_binding
      real(kind=c_double) :: timestepLengthLimit
    end subroutine precicef_advance

    subroutine precicef_finalize() bind(c, name='precicef_finalize_')

      use, intrinsic :: iso_c_binding
    end subroutine precicef_finalize

    subroutine precicef_get_dims(dimensions) &
      &  bind(c, name='precicef_get_dims_')

      use, intrinsic :: iso_c_binding
      integer(kind=c_int) :: dimensions
    end subroutine precicef_get_dims

    subroutine precicef_ongoing(isOngoing) &
      &  bind(c, name='precicef_ongoing_')

      use, intrinsic :: iso_c_binding
      integer(kind=c_int) :: isOngoing
    end subroutine precicef_ongoing

    subroutine precicef_write_data_required(computedTimestepLength, &
      &                                     isRequired) &
      &  bind(c, name='precicef_write_data_required_')

      use, intrinsic :: iso_c_binding
      real(kind=c_double) :: computedTimestepLength
      integer(kind=c_int) :: isRequired
    end subroutine precicef_write_data_required

    subroutine precicef_read_data_available(isAvailable) &
      &  bind(c, name='precicef_read_data_available_')

      use, intrinsic :: iso_c_binding
      integer(kind=c_int) :: isAvailable
    end subroutine precicef_read_data_available

    subroutine precicef_action_required(action, isRequired, lengthAction) &
      &  bind(c, name='precicef_action_required_')

      use, intrinsic :: iso_c_binding
      character(kind=c_char), dimension(*) :: action
      integer(kind=c_int)                  :: isRequired
      integer(kind=c_int), value           :: lengthAction
    end subroutine precicef_action_required

    subroutine precicef_fulfilled_action(action, lengthAction) &
      &  bind(c, name='precicef_fulfilled_action_')

      use, intrinsic :: iso_c_binding
      character(kind=c_char), dimension(*) :: action
      integer(kind=c_int), value           :: lengthAction
    end subroutine precicef_fulfilled_action

    subroutine precicef_get_mesh_id(meshName, meshID, lengthMeshName) &
      &  bind(c, name='precicef_get_mesh_id_')

      use, intrinsic :: iso_c_binding
      character(kind=c_char), dimension(*) :: meshName
      integer(kind=c_int)                  :: meshID
      integer(kind=c_int), value           :: lengthMeshName
    end subroutine precicef_get_mesh_id

    subroutine precicef_has_data(dataName, hasData, lengthDataName, meshID) &
      &  bind(c, name='precicef_has_data_')

      use, intrinsic :: iso_c_binding
      character(kind=c_char), dimension(*) :: dataName
      integer(kind=c_int)                  :: hasData
      integer(kind=c_int), value           :: lengthDataName
      integer(kind=c_int)                  :: meshID
    end subroutine precicef_has_data

    subroutine precicef_get_data_id(dataName, meshID,dataID, lengthDataName ) &
      &  bind(c, name='precicef_get_data_id_')

      use, intrinsic :: iso_c_binding
      character(kind=c_char), dimension(*) :: dataName
      integer(kind=c_int)                  :: meshID
      integer(kind=c_int)                  :: dataID
      integer(kind=c_int), value           :: lengthDataName
    end subroutine precicef_get_data_id

    subroutine precicef_set_vertices(meshID, meshsize, positions, positionIDs) &
      &  bind(c, name='precicef_set_vertices_')

      use, intrinsic :: iso_c_binding
      integer(kind=c_int) :: meshID
      integer(kind=c_int) :: meshsize
      real(kind=c_double) :: positions(*)
      integer(kind=c_int) :: positionIDs(*)
    end subroutine precicef_set_vertices

    subroutine precicef_set_vertex(meshID, position, vertexID) &
      &  bind(c, name='precicef_set_vertex_')

      use, intrinsic :: iso_c_binding
      integer(kind=c_int) :: meshID
      real(kind=c_double) :: position(3)
      integer(kind=c_int) :: vertexID
    end subroutine precicef_set_vertex

    subroutine precicef_set_edge(meshID, firstVertexID, secondVertexID, &
      &                          edgeID) &
      &  bind(c, name='precicef_set_edge_')

      use, intrinsic :: iso_c_binding
      integer(kind=c_int) :: meshID
      integer(kind=c_int) :: firstVertexID
      integer(kind=c_int) :: secondVertexID
      integer(kind=c_int) :: edgeID
    end subroutine precicef_set_edge

    subroutine precicef_set_triangle(meshID, firstEdgeID, secondEdgeID, &
      &                              thirdEdgeID) &
      &  bind(c, name='precicef_set_triangle_')

      use, intrinsic :: iso_c_binding
      integer(kind=c_int) :: meshID
      integer(kind=c_int) :: firstEdgeID
      integer(kind=c_int) :: secondEdgeID
      integer(kind=c_int) :: thirdEdgeID
    end subroutine precicef_set_triangle
 
    subroutine precicef_read_sdata( dataID, valueIndex, dataValue) &
      &  bind(c, name='precicef_read_sdata_')

      use, intrinsic :: iso_c_binding
      integer(kind=c_int) :: dataID
      integer(kind=c_int) :: valueIndex
      real(kind=c_double) :: dataValue   
    end subroutine precicef_read_sdata

    subroutine precicef_read_bsdata( dataID, blocksize, valueIndices, values) &
      &  bind(c, name='precicef_read_bsdata_')

      use, intrinsic :: iso_c_binding
      integer(kind=c_int) :: dataID
      integer(kind=c_int) :: blocksize
      integer(kind=c_int) :: valueIndices(*)
      real(kind=c_double) :: values(*)   
    end subroutine precicef_read_bsdata

    subroutine precicef_read_vdata( dataID, valueIndex, dataValue) &
      &  bind(c, name='precicef_read_vdata_')

      use, intrinsic :: iso_c_binding
      integer(kind=c_int) :: dataID
      integer(kind=c_int) :: valueIndex
      real(kind=c_double) :: dataValue(*)
    end subroutine precicef_read_vdata

    subroutine precicef_write_sdata( dataID, valueIndex, dataValue) &
      &  bind(c, name='precicef_write_sdata_')

      use, intrinsic :: iso_c_binding
      integer(kind=c_int) :: dataID
      integer(kind=c_int) :: valueIndex
      real(kind=c_double) :: dataValue   
    end subroutine precicef_write_sdata

    subroutine precicef_write_bsdata( dataID, blockSize, valueIndices, values) &
      &  bind(c, name='precicef_write_bsdata_')

      use, intrinsic :: iso_c_binding
      integer(kind=c_int) :: dataID
      integer(kind=c_int) :: blockSize
      integer(kind=c_int) :: valueIndices(*)
      real(kind=c_double) :: values(*)   
    end subroutine precicef_write_bsdata

    subroutine precicef_write_vdata( dataID, valueIndex, dataValue) &
      &  bind(c, name='precicef_write_vdata_')

      use, intrinsic :: iso_c_binding
      integer(kind=c_int) :: dataID
      integer(kind=c_int) :: valueIndex
      real(kind=c_double), dimension(*) :: dataValue 
    end subroutine precicef_write_vdata

    subroutine precicef_write_data_available(isAvailable) &
      &  bind(c, name='precicef_write_data_available_')

      use, intrinsic :: iso_c_binding
      integer(kind=c_int) :: isAvailable
    end subroutine precicef_write_data_available
    
    subroutine precicef_initialize_data() &
      &  bind(c, name='precicef_initialize_data_')

      use, intrinsic :: iso_c_binding
    end subroutine precicef_initialize_data

    !!! TO BE CONTINUED ...

  end interface

end module PreCICE_solver_if_module
