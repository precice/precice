subroutine precice_main_f2c_nsd_create_client_instance(self,host,port,buffer_size) bind(c)
     use, intrinsic :: iso_c_binding
     integer(kind=c_long_long)::self
     type(c_ptr)::host
     integer(kind=c_int)::port
     integer(kind=c_int)::buffer_size
end subroutine precice_main_f2c_nsd_create_client_instance


subroutine precice_main_f2c_nsd_destroy_instance(self) bind(c)
     use, intrinsic :: iso_c_binding
     integer(kind=c_long_long)::self
end subroutine precice_main_f2c_nsd_destroy_instance

subroutine  precice_main_f2c_nsd_main(self) bind(c)
     use, intrinsic :: iso_c_binding
     integer(kind=c_long_long)::self
     
end subroutine precice_main_f2c_nsd_main
