subroutine  fsi_fsicommc2socket_plain_port_create_client_instance(self,host,port,buffer_size) bind(c)
     use, intrinsic :: iso_c_binding
     integer(kind=c_long_long)::self
     type(c_ptr)::host
     integer(kind=c_int)::port
     integer(kind=c_int)::buffer_size
     
     
end subroutine fsi_fsicommc2socket_plain_port_create_client_instance


subroutine  fsi_fsicommc2socket_plain_port_destroy_instance(self) bind(c)
     use, intrinsic :: iso_c_binding
     integer(kind=c_long_long)::self
end subroutine fsi_fsicommc2socket_plain_port_destroy_instance

subroutine  fsi_fsicommc2socket_plain_port_transferCoordinates(self,&
	coordId,coordId_len,&
	offsets,offsets_len,&
	hosts,hosts_len) bind(c)
     use, intrinsic :: iso_c_binding
     integer(kind=c_long_long)::self
     	integer(kind=c_int),intent(in),dimension(*)::coordId
	integer(kind=c_int),intent(in)::coordId_len
	integer(kind=c_int),intent(in),dimension(*)::offsets
	integer(kind=c_int),intent(in)::offsets_len
	integer(kind=c_int),intent(in)::hosts_len
	type(c_ptr),dimension(*),intent(in)::hosts

end subroutine fsi_fsicommc2socket_plain_port_transferCoordinates
subroutine  fsi_fsicommc2socket_plain_port_startDataTransfer(self) bind(c)
     use, intrinsic :: iso_c_binding
     integer(kind=c_long_long)::self
     
end subroutine fsi_fsicommc2socket_plain_port_startDataTransfer
subroutine  fsi_fsicommc2socket_plain_port_endDataTransfer(self,&
	ack) bind(c)
     use, intrinsic :: iso_c_binding
     integer(kind=c_long_long)::self
     	integer(kind=c_int),intent(inout)::ack

end subroutine fsi_fsicommc2socket_plain_port_endDataTransfer
