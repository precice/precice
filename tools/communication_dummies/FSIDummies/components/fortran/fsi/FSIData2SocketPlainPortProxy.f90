subroutine  fsi_fsidatac2socket_plain_port_create_client_instance(self,host,port,buffer_size) bind(c)
     use, intrinsic :: iso_c_binding
     integer(kind=c_long_long)::self
     type(c_ptr)::host
     integer(kind=c_int)::port
     integer(kind=c_int)::buffer_size
     
     
end subroutine fsi_fsidatac2socket_plain_port_create_client_instance


subroutine  fsi_fsidatac2socket_plain_port_destroy_instance(self) bind(c)
     use, intrinsic :: iso_c_binding
     integer(kind=c_long_long)::self
end subroutine fsi_fsidatac2socket_plain_port_destroy_instance

subroutine  fsi_fsidatac2socket_plain_port_transferData(self,&
	coordId,coordId_len,&
	data,data_len) bind(c)
     use, intrinsic :: iso_c_binding
     integer(kind=c_long_long)::self
     	integer(kind=c_int),intent(in),dimension(*)::coordId
	integer(kind=c_int),intent(in)::coordId_len
	real(kind=c_double),intent(in),dimension(*)::data
	integer(kind=c_int),intent(in)::data_len

end subroutine fsi_fsidatac2socket_plain_port_transferData
subroutine  fsi_fsidatac2socket_plain_port_dataAck(self,&
	ack) bind(c)
     use, intrinsic :: iso_c_binding
     integer(kind=c_long_long)::self
     	integer(kind=c_int),intent(inout)::ack

end subroutine fsi_fsidatac2socket_plain_port_dataAck
