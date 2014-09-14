subroutine  precice_initializerc2socket_plain_port_create_client_instance(self,host,port,buffer_size) bind(c)
     use, intrinsic :: iso_c_binding
     integer(kind=c_long_long)::self
     type(c_ptr)::host
     integer(kind=c_int)::port
     integer(kind=c_int)::buffer_size
     
     
end subroutine precice_initializerc2socket_plain_port_create_client_instance


subroutine  precice_initializerc2socket_plain_port_destroy_instance(self) bind(c)
     use, intrinsic :: iso_c_binding
     integer(kind=c_long_long)::self
end subroutine precice_initializerc2socket_plain_port_destroy_instance

subroutine  precice_initializerc2socket_plain_port_initializeAddresses(self,&
	addresses,addresses_len) bind(c)
     use, intrinsic :: iso_c_binding
     integer(kind=c_long_long)::self
     	integer(kind=c_int),intent(in)::addresses_len
	type(c_ptr),dimension(*),intent(in)::addresses

end subroutine precice_initializerc2socket_plain_port_initializeAddresses
subroutine  precice_initializerc2socket_plain_port_initializeVertexes(self,&
	vertexes,vertexes_len) bind(c)
     use, intrinsic :: iso_c_binding
     integer(kind=c_long_long)::self
     	integer(kind=c_int),intent(in),dimension(*)::vertexes
	integer(kind=c_int),intent(in)::vertexes_len

end subroutine precice_initializerc2socket_plain_port_initializeVertexes
