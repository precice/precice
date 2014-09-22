module precice_Communicator2SocketPort 
use, intrinsic :: iso_c_binding
implicit none


type, public :: Communicator2SocketPort
     
     integer(kind=c_long_long )::reference
     contains
     procedure,public::create_port_client_instance
     procedure,public::create_port_client_instance_for_c
     procedure,public::destroyPortInstance
     procedure,public::setData

end type Communicator2SocketPort
contains

subroutine create_port_client_instance_for_c(this,host,port,buffer_size)
    class(Communicator2SocketPort)::this
    integer(kind=c_int)::port
    character(kind=c_char),dimension(*)::host
    integer(kind=c_int)::buffer_size
    call precice_communicatorc2socket_plain_port_create_client_instance(this%reference,host,port,buffer_size)
    

end subroutine create_port_client_instance_for_c

subroutine create_port_client_instance(this,host,port,buffer_size)
    class(Communicator2SocketPort)::this
    integer::port
    character(*)::host
    integer::buffer_size
    call this%create_port_client_instance_for_c(host//c_null_char,port,buffer_size)
    

end subroutine create_port_client_instance



subroutine destroyPortInstance(this)
     class(Communicator2SocketPort)::this
     call precice_communicatorc2socket_plain_port_destroy_instance(this%reference)

end subroutine destroyPortInstance

subroutine setData(this,&
	data,&
	index,&
	rank,&
	tag)
     use, intrinsic :: iso_c_binding
     class(Communicator2SocketPort)::this
     	real(8),intent(in)::data
	integer,intent(in)::index
	integer,intent(in)::rank
	integer,intent(inout)::tag

     
     call precice_communicatorc2socket_plain_port_setData(this%reference,&
data,&
index,&
rank,&
tag)
end subroutine setData

end module  precice_Communicator2SocketPort
