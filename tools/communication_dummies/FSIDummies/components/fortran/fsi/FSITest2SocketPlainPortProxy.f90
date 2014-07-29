subroutine  fsi_fsitestc2socket_plain_port_create_client_instance(self,host,port,buffer_size) bind(c)
     use, intrinsic :: iso_c_binding
     integer(kind=c_long_long)::self
     type(c_ptr)::host
     integer(kind=c_int)::port
     integer(kind=c_int)::buffer_size
     
     
end subroutine fsi_fsitestc2socket_plain_port_create_client_instance


subroutine  fsi_fsitestc2socket_plain_port_destroy_instance(self) bind(c)
     use, intrinsic :: iso_c_binding
     integer(kind=c_long_long)::self
end subroutine fsi_fsitestc2socket_plain_port_destroy_instance

subroutine  fsi_fsitestc2socket_plain_port_test(self) bind(c)
     use, intrinsic :: iso_c_binding
     integer(kind=c_long_long)::self
     
end subroutine fsi_fsitestc2socket_plain_port_test
