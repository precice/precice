subroutine precice_initializer_f2c_nsd_create_client_instance(self,host,port,buffer_size) bind(c)
     use, intrinsic :: iso_c_binding
     integer(kind=c_long_long)::self
     type(c_ptr)::host
     integer(kind=c_int)::port
     integer(kind=c_int)::buffer_size
end subroutine precice_initializer_f2c_nsd_create_client_instance


subroutine precice_initializer_f2c_nsd_destroy_instance(self) bind(c)
     use, intrinsic :: iso_c_binding
     integer(kind=c_long_long)::self
end subroutine precice_initializer_f2c_nsd_destroy_instance

subroutine  precice_initializer_f2c_nsd_acknowledge(self,&
	identifier,&
	tag) bind(c)
     use, intrinsic :: iso_c_binding
     integer(kind=c_long_long)::self
     	integer(kind=c_int),intent(in)::identifier
	integer(kind=c_int),intent(inout)::tag

end subroutine precice_initializer_f2c_nsd_acknowledge
subroutine  precice_initializer_f2c_nsd_initialize(self,&
	addresses,addresses_len,&
	vertexes,vertexes_len) bind(c)
     use, intrinsic :: iso_c_binding
     integer(kind=c_long_long)::self
     	integer(kind=c_int),intent(in)::addresses_len
	type(c_ptr),dimension(*),intent(in)::addresses
	integer(kind=c_int),intent(in),dimension(*)::vertexes
	integer(kind=c_int),intent(in)::vertexes_len

end subroutine precice_initializer_f2c_nsd_initialize
