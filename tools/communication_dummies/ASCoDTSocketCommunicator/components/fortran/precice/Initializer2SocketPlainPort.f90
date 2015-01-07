module precice_Initializer2SocketPort 
use, intrinsic :: iso_c_binding
implicit none


type, public :: Initializer2SocketPort
     
     integer(kind=c_long_long )::reference
     contains
     procedure,public::create_port_client_instance
     procedure,public::create_port_client_instance_for_c
     procedure,public::destroyPortInstance
     procedure,public::acknowledge
procedure,public::initialize

end type Initializer2SocketPort
contains

subroutine create_port_client_instance_for_c(this,host,port,buffer_size)
    class(Initializer2SocketPort)::this
    integer(kind=c_int)::port
    character(kind=c_char),dimension(*)::host
    integer(kind=c_int)::buffer_size
    call precice_initializerc2socket_plain_port_create_client_instance(this%reference,host,port,buffer_size)
    

end subroutine create_port_client_instance_for_c

subroutine create_port_client_instance(this,host,port,buffer_size)
    class(Initializer2SocketPort)::this
    integer::port
    character(*)::host
    integer::buffer_size
    call this%create_port_client_instance_for_c(host//c_null_char,port,buffer_size)
    

end subroutine create_port_client_instance



subroutine destroyPortInstance(this)
     class(Initializer2SocketPort)::this
     call precice_initializerc2socket_plain_port_destroy_instance(this%reference)

end subroutine destroyPortInstance

subroutine initialize(this,&
	addresses,addresses_len,&
	vertexes,vertexes_len)
     use, intrinsic :: iso_c_binding
     class(Initializer2SocketPort)::this
     	character(*),intent(in),dimension(*)::addresses
	integer,intent(in)::addresses_len
	type(c_ptr),dimension(addresses_len) :: addressesPtrArray
	integer::addresses_ns
	character(255), dimension(addresses_len), target :: addressesFSArray
	integer,intent(in),dimension(*)::vertexes
	integer,intent(in)::vertexes_len
	do addresses_ns = 1, addresses_len
		addressesFSArray(addresses_ns) = addresses(addresses_ns)// C_NULL_CHAR
		addressesPtrArray(addresses_ns) = C_LOC(addressesFSArray(addresses_ns))
	end do

     
     call precice_initializerc2socket_plain_port_initialize(this%reference,&
addresses,addresses_len,&
vertexes,vertexes_len)
end subroutine initialize
subroutine acknowledge(this,&
	identifier,&
	tag)
     use, intrinsic :: iso_c_binding
     class(Initializer2SocketPort)::this
     	integer,intent(in)::identifier
	integer,intent(inout)::tag

     
     call precice_initializerc2socket_plain_port_acknowledge(this%reference,&
identifier,&
tag)
end subroutine acknowledge

end module  precice_Initializer2SocketPort
