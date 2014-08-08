module fsi_FSIData2SocketPort 
use, intrinsic :: iso_c_binding
implicit none


type, public :: FSIData2SocketPort
     
     integer(kind=c_long_long )::reference
     contains
     procedure,public::create_port_client_instance
     procedure,public::create_port_client_instance_for_c
     procedure,public::destroyPortInstance
     procedure,public::transferData
procedure,public::dataAck

end type FSIData2SocketPort
contains

subroutine create_port_client_instance_for_c(this,host,port,buffer_size)
    class(FSIData2SocketPort)::this
    integer(kind=c_int)::port
    character(kind=c_char),dimension(*)::host
    integer(kind=c_int)::buffer_size
    call fsi_fsidatac2socket_plain_port_create_client_instance(this%reference,host,port,buffer_size)
    

end subroutine create_port_client_instance_for_c

subroutine create_port_client_instance(this,host,port,buffer_size)
    class(FSIData2SocketPort)::this
    integer::port
    character(*)::host
    integer::buffer_size
    call this%create_port_client_instance_for_c(host//c_null_char,port,buffer_size)
    

end subroutine create_port_client_instance



subroutine destroyPortInstance(this)
     class(FSIData2SocketPort)::this
     call fsi_fsidatac2socket_plain_port_destroy_instance(this%reference)

end subroutine destroyPortInstance

subroutine dataAck(this,&
	ack)
     use, intrinsic :: iso_c_binding
     class(FSIData2SocketPort)::this
     	integer,intent(inout)::ack

     
     call fsi_fsidatac2socket_plain_port_dataAck(this%reference,&
ack)
end subroutine dataAck
subroutine transferData(this,&
	coordId,coordId_len,&
	data,data_len)
     use, intrinsic :: iso_c_binding
     class(FSIData2SocketPort)::this
     	integer,intent(in),dimension(*)::coordId
	integer,intent(in)::coordId_len
	real(8),intent(in),dimension(*)::data
	integer,intent(in)::data_len

     
     call fsi_fsidatac2socket_plain_port_transferData(this%reference,&
coordId,coordId_len,&
data,data_len)
end subroutine transferData

end module  fsi_FSIData2SocketPort
