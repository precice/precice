module fsi_FSIDataFNativeSocketDispatcher
use, intrinsic :: iso_c_binding
implicit none


type, public :: FSIDataNativeSocketDispatcher
     integer(kind=c_long_long )::reference
     contains
     procedure,public::createClientDispatcherInstanceForC
     procedure,public::createClientDispatcherInstance
     procedure,public::destroyDispatcherInstance
     
     
     	procedure,public::transferData
	procedure,private::transferData_internal
	procedure,public::dataAck
	procedure,private::dataAck_internal

end type FSIDataNativeSocketDispatcher
contains
subroutine createClientDispatcherInstanceForC(this,host,port,buffer_size)
    class(FSIDataNativeSocketDispatcher)::this
    character(kind=c_char),dimension(*)::host
    integer(kind=c_int)::port
    integer(kind=c_int)::buffer_size
    this%reference=0
    call fsi_fsidata_f2c_nsd_create_client_instance(this%reference,host,port,buffer_size)
    
    

end subroutine createClientDispatcherInstanceForC

subroutine createClientDispatcherInstance(this,host,port,buffer_size)
    class(FSIDataNativeSocketDispatcher)::this
    character(*)::host
    integer::port
    integer::buffer_size
    call this%createClientDispatcherInstanceForC(host//c_null_char,port,buffer_size)
    
    

end subroutine createClientDispatcherInstance

subroutine destroyDispatcherInstance(this)
     class(FSIDataNativeSocketDispatcher)::this
     call fsi_fsidata_f2c_nsd_destroy_instance(this%reference)

end subroutine destroyDispatcherInstance

subroutine dataAck_internal(this,&
	ack)
     use, intrinsic :: iso_c_binding
     class(FSIDataNativeSocketDispatcher)::this
     	integer(kind=c_int),intent(inout)::ack

     call fsi_fsidata_f2c_nsd_dataAck(this%reference,&
ack)
end subroutine dataAck_internal

subroutine dataAck(this,&
	ack)
     use, intrinsic :: iso_c_binding
     class(FSIDataNativeSocketDispatcher)::this
     	integer,intent(inout)::ack

     call this%dataAck_internal(&
ack)
end subroutine dataAck
subroutine transferData_internal(this,&
	coordId,coordId_len,&
	data,data_len)
     use, intrinsic :: iso_c_binding
     class(FSIDataNativeSocketDispatcher)::this
     	integer(kind=c_int),intent(in),dimension(*)::coordId
	integer(kind=c_int),intent(in)::coordId_len
	real(kind=c_double),intent(in),dimension(*)::data
	integer(kind=c_int),intent(in)::data_len

     call fsi_fsidata_f2c_nsd_transferData(this%reference,&
coordId,coordId_len,&
data,data_len)
end subroutine transferData_internal

subroutine transferData(this,&
	coordId,coordId_len,&
	data,data_len)
     use, intrinsic :: iso_c_binding
     class(FSIDataNativeSocketDispatcher)::this
     	integer,intent(in),dimension(*)::coordId
	integer,intent(in)::coordId_len
	real(8),intent(in),dimension(*)::data
	integer,intent(in)::data_len

     call this%transferData_internal(&
coordId,coordId_len,&
data,data_len)
end subroutine transferData

end module  fsi_FSIDataFNativeSocketDispatcher
