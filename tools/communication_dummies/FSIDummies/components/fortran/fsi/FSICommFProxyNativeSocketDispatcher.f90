subroutine fsi_fsicomm_f2c_nsd_create_client_instance(self,host,port,buffer_size) bind(c)
     use, intrinsic :: iso_c_binding
     integer(kind=c_long_long)::self
     type(c_ptr)::host
     integer(kind=c_int)::port
     integer(kind=c_int)::buffer_size
end subroutine fsi_fsicomm_f2c_nsd_create_client_instance


subroutine fsi_fsicomm_f2c_nsd_destroy_instance(self) bind(c)
     use, intrinsic :: iso_c_binding
     integer(kind=c_long_long)::self
end subroutine fsi_fsicomm_f2c_nsd_destroy_instance

subroutine  fsi_fsicomm_f2c_nsd_transferCoordinates(self,&
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

end subroutine fsi_fsicomm_f2c_nsd_transferCoordinates
subroutine  fsi_fsicomm_f2c_nsd_transferData(self,&
	data,data_len) bind(c)
     use, intrinsic :: iso_c_binding
     integer(kind=c_long_long)::self
     	real(kind=c_double),intent(in),dimension(*)::data
	integer(kind=c_int),intent(in)::data_len

end subroutine fsi_fsicomm_f2c_nsd_transferData
