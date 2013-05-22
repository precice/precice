module PreCICE_solver_if_module
  use, intrinsic :: iso_c_binding
  implicit none

  interface

    subroutine precicef_create(accessorName, configFileName, &
      &                        solverProcessIndex, solverProcessSize, &
      &                        lengthAccessorName, lengthConfigFileName) &
      &  bind(c, name='precicef_create_')

      use, intrinsic :: iso_c_binding
      character(kind=c_char), dimension(*) :: accessorName
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

!TODO!    subroutine precicef_get_dims(dimensions) &
!TODO!      &  bind(c, name='precicef_get_dims_')
!TODO!
!TODO!      use, intrinsic :: iso_c_binding
!TODO!      integer(kind=c_int) :: dimensions
!TODO!    end subroutine precicef_get_dims
!TODO!
!TODO!    subroutine precicef_ongoing(isOngoing) &
!TODO!      &  bind(c, name='precicef_ongoing_')
!TODO!
!TODO!      use, intrinsic :: iso_c_binding
!TODO!      integer(kind=c_int) :: isOngoing
!TODO!    end subroutine precicef_ongoing
!TODO!
!TODO!    subroutine precicef_write_data_required(computedTimestepLength, &
!TODO!      &                                     isRequired) &
!TODO!      &  bind(c, name='precicef_write_data_required_')
!TODO!
!TODO!      use, intrinsic :: iso_c_binding
!TODO!      real(kind=c_double) :: computedTimestepLength
!TODO!      integer(kind=c_int) :: isRequired
!TODO!    end subroutine precicef_write_data_required
!TODO!
!TODO!    subroutine precicef_read_data_available(isAvailable) &
!TODO!      &  bind(c, name='precicef_read_data_available_')
!TODO!
!TODO!      use, intrinsic :: iso_c_binding
!TODO!      integer(kind=c_int) :: isAvailable
!TODO!    end subroutine precicef_read_data_available
!TODO!
!TODO!    subroutine precicef_action_required(action, isRequired, lengthAction) &
!TODO!      &  bind(c, name='precicef_action_required_')
!TODO!
!TODO!      use, intrinsic :: iso_c_binding
!TODO!      character(kind=c_char), dimension(*) :: action
!TODO!      integer(kind=c_int) :: isRequired
!TODO!      integer(kind=c_int), value :: lengthAction
!TODO!    end subroutine precicef_action_required
!TODO!
!TODO!    subroutine precicef_fulfilled_action(action, lengthAction) &
!TODO!      &  bind(c, name='precicef_fulfilled_action_')
!TODO!
!TODO!      use, intrinsic :: iso_c_binding
!TODO!      character(kind=c_char), dimension(*) :: action
!TODO!      integer(kind=c_int), value :: lengthAction
!TODO!    end subroutine precicef_fulfilled_action
!TODO!
!TODO!    subroutine precicef_get_mesh_id(geometryName, meshID, lengthGeometryName) &
!TODO!      &  bind(c, name='precicef_get_mesh_id_')
!TODO!
!TODO!      use, intrinsic :: iso_c_binding
!TODO!      character(kind=c_char), dimension(*) :: geometryName
!TODO!      integer(kind=c_int) :: meshID
!TODO!      integer(kind=c_int), value :: lengthGeometryName
!TODO!    end subroutine precicef_get_mesh_id
!TODO!
!TODO!    subroutine precicef_has_data(dataName, hasData, lengthDataName) &
!TODO!      &  bind(c, name='precicef_has_data_')
!TODO!
!TODO!      use, intrinsic :: iso_c_binding
!TODO!      character(kind=c_char), dimension(*) :: dataName
!TODO!      integer(kind=c_int) :: hasData
!TODO!      integer(kind=c_int), value :: lengthDataName
!TODO!    end subroutine precicef_has_data
!TODO!
!TODO!    subroutine precicef_get_data_id(dataName, dataID, lengthDataName) &
!TODO!      &  bind(c, name='precicef_get_data_id_')
!TODO!
!TODO!      use, intrinsic :: iso_c_binding
!TODO!      character(kind=c_char), dimension(*) :: dataName
!TODO!      integer(kind=c_int) :: dataID
!TODO!      integer(kind=c_int), value :: lengthDataName
!TODO!    end subroutine precicef_get_data_id
!TODO!
!TODO!    subroutine precicef_set_vertex(meshID, position, vertexID) &
!TODO!      &  bind(c, name='precicef_set_vertex_')
!TODO!
!TODO!      use, intrinsic :: iso_c_binding
!TODO!      integer(kind=c_int) :: meshID
!TODO!      real(kind=c_double) :: position
!TODO!      integer(kind=c_int) :: vertexID
!TODO!    end subroutine precicef_set_vertex
!TODO!
!TODO!    subroutine precicef_set_read_pos(meshID, position, vertexID) &
!TODO!      &  bind(c, name='precicef_set_read_pos_')
!TODO!
!TODO!      use, intrinsic :: iso_c_binding
!TODO!      integer(kind=c_int) :: meshID
!TODO!      real(kind=c_double) :: position
!TODO!      integer(kind=c_int) :: positionID
!TODO!    end subroutine precicef_set_read_pos
!TODO!
!TODO!    subroutine precicef_write_pos(meshID, position, vertexID) &
!TODO!      &  bind(c, name='precicef_set_write_pos_')
!TODO!
!TODO!      use, intrinsic :: iso_c_binding
!TODO!      integer(kind=c_int) :: meshID
!TODO!      real(kind=c_double) :: position
!TODO!      integer(kind=c_int) :: positionID
!TODO!    end subroutine precicef_set_write_pos
!TODO!
!TODO!    subroutine precicef_set_edge(meshID, firstVertexID, secondVertexID, &
!TODO!      &                          edgeID) &
!TODO!      &  bind(c, name='precicef_set_edge_')
!TODO!
!TODO!      use, intrinsic :: iso_c_binding
!TODO!      integer(kind=c_int) :: meshID
!TODO!      integer(kind=c_int) :: firstVertexID
!TODO!      integer(kind=c_int) :: secondVertexID
!TODO!      integer(kind=c_int) :: edgeID
!TODO!    end subroutine precicef_set_edge
!TODO!
!TODO!    subroutine precicef_set_triangle(meshID, firstEdgeID, secondEdgeID, &
!TODO!      &                              thirdEdgeID) &
!TODO!      &  bind(c, name='precicef_set_edge_')
!TODO!
!TODO!      use, intrinsic :: iso_c_binding
!TODO!      integer(kind=c_int) :: meshID
!TODO!      integer(kind=c_int) :: firstEdgeID
!TODO!      integer(kind=c_int) :: secondEdgeID
!TODO!      integer(kind=c_int) :: thirdEdgeID
!TODO!    end subroutine precicef_set_triangle

    !!! TO BE CONTINUED ...

  end interface

end module PreCICE_solver_if_module
