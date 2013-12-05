module PreCICE_consts_module
  use, intrinsic :: iso_c_binding

  implicit none

  interface

    subroutine precicef_name_config(nameConfig, lengthNameConfig) &
      &        bind(c, name="precicef_name_config_")
      use, intrinsic :: iso_c_binding
      character(kind=c_char), dimension(*) :: nameConfig
      integer(kind=c_int), value :: lengthNameConfig
    end subroutine precicef_name_config


    subroutine precicef_action_write_iter_checkp(nameAction, lengthNameAction) &
      &        bind(c, name="precicef_action_write_iter_checkp_")
      use, intrinsic :: iso_c_binding
      character(kind=c_char), dimension(*) :: nameAction
      integer(kind=c_int), value :: lengthNameAction
    end subroutine precicef_action_write_iter_checkp

    subroutine precicef_action_write_initial_data(nameAction, lengthNameAction) &
      &        bind(c, name="precicef_action_write_initial_data_")
      use, intrinsic :: iso_c_binding
      character(kind=c_char), dimension(*) :: nameAction
      integer(kind=c_int), value :: lengthNameAction
    end subroutine precicef_action_write_initial_data

    subroutine precicef_action_read_iter_checkp(nameAction, lengthNameAction) &
      &        bind(c, name="precicef_action_read_iter_checkp_")
      use, intrinsic :: iso_c_binding
      character(kind=c_char), dimension(*) :: nameAction
      integer(kind=c_int), value :: lengthNameAction
    end subroutine precicef_action_read_iter_checkp


    subroutine precicef_action_write_sim_checkp(nameAction, lengthNameAction) &
      &        bind(c, name="precicef_action_write_sim_checkp_")
      use, intrinsic :: iso_c_binding
      character(kind=c_char), dimension(*) :: nameAction
      integer(kind=c_int), value :: lengthNameAction
    end subroutine precicef_action_write_sim_checkp

    subroutine precicef_action_read_sim_checkp(nameAction, lengthNameAction) &
      &        bind(c, name="precicef_action_read_sim_checkp_")
      use, intrinsic :: iso_c_binding
      character(kind=c_char), dimension(*) :: nameAction
      integer(kind=c_int), value :: lengthNameAction
    end subroutine precicef_action_read_sim_checkp

  end interface

end module PreCICE_consts_module
