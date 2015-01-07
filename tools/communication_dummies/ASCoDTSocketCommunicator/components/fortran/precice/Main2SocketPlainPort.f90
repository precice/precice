module precice_Main2SocketPort 
use, intrinsic :: iso_c_binding
implicit none


type, public :: Main2SocketPort
     
     integer(kind=c_long_long )::reference
     contains
     procedure,public::create_port_client_instance
     procedure,public::create_port_client_instance_for_c
     procedure,public::destroyPortInstance
     procedure,public::main

end type Main2SocketPort
contains

subroutine create_port_client_instance_for_c(this,host,port,buffer_size)
    class(Main2SocketPort)::this
    integer(kind=c_int)::port
    character(kind=c_char),dimension(*)::host
    integer(kind=c_int)::buffer_size
    call precice_mainc2socket_plain_port_create_client_instance(this%reference,host,port,buffer_size)
    

end subroutine create_port_client_instance_for_c

subroutine create_port_client_instance(this,host,port,buffer_size)
    class(Main2SocketPort)::this
    integer::port
    character(*)::host
    integer::buffer_size
    call this%create_port_client_instance_for_c(host//c_null_char,port,buffer_size)
    

end subroutine create_port_client_instance



subroutine destroyPortInstance(this)
     class(Main2SocketPort)::this
     call precice_mainc2socket_plain_port_destroy_instance(this%reference)

end subroutine destroyPortInstance

subroutine main(this)
     use, intrinsic :: iso_c_binding
     class(Main2SocketPort)::this
     
     
     call precice_mainc2socket_plain_port_main(this%reference)
end subroutine main

end module  precice_Main2SocketPort
