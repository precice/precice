module fsi_FSIComm2SocketPort 
use, intrinsic :: iso_c_binding
implicit none


type, public :: FSIComm2SocketPort
     
     integer(kind=c_long_long )::reference
     contains
     procedure,public::create_port_client_instance
     procedure,public::create_port_client_instance_for_c
     procedure,public::destroyPortInstance
     procedure,public::transferCoordinates
procedure,public::transferData

end type FSIComm2SocketPort
contains

subroutine create_port_client_instance_for_c(this,host,port,buffer_size)
    class(FSIComm2SocketPort)::this
    integer(kind=c_int)::port
    character(kind=c_char),dimension(*)::host
    integer(kind=c_int)::buffer_size
    call fsi_fsicommc2socket_plain_port_create_client_instance(this%reference,host,port,buffer_size)
    

end subroutine create_port_client_instance_for_c

subroutine create_port_client_instance(this,host,port,buffer_size)
    class(FSIComm2SocketPort)::this
    integer::port
    character(*)::host
    integer::buffer_size
    call this%create_port_client_instance_for_c(host//c_null_char,port,buffer_size)
    

end subroutine create_port_client_instance



subroutine destroyPortInstance(this)
     class(FSIComm2SocketPort)::this
     call fsi_fsicommc2socket_plain_port_destroy_instance(this%reference)

end subroutine destroyPortInstance

subroutine transferData(this,&
	data,data_len)
     use, intrinsic :: iso_c_binding
     class(FSIComm2SocketPort)::this
     	real(8),intent(in),dimension(*)::data
	integer,intent(in)::data_len

     
     call fsi_fsicommc2socket_plain_port_transferData(this%reference,&
data,data_len)
end subroutine transferData
subroutine transferCoordinates(this,&
	coordId,coordId_len,&
	offsets,offsets_len,&
	hosts,hosts_len)
     use, intrinsic :: iso_c_binding
     class(FSIComm2SocketPort)::this
     	integer,intent(in),dimension(*)::coordId
	integer,intent(in)::coordId_len
	integer,intent(in),dimension(*)::offsets
	integer,intent(in)::offsets_len
	character(*),intent(in),dimension(*)::hosts
	integer,intent(in)::hosts_len
	type(c_ptr),dimension(hosts_len) :: hostsPtrArray
	integer::hosts_ns
	character(255), dimension(hosts_len), target :: hostsFSArray
	do hosts_ns = 1, hosts_len
		hostsFSArray(hosts_ns) = hosts(hosts_ns)// C_NULL_CHAR
		hostsPtrArray(hosts_ns) = C_LOC(hostsFSArray(hosts_ns))
	end do

     
     call fsi_fsicommc2socket_plain_port_transferCoordinates(this%reference,&
coordId,coordId_len,&
offsets,offsets_len,&
hosts,hosts_len)
end subroutine transferCoordinates

end module  fsi_FSIComm2SocketPort
