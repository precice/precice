module fsi_FSICommFNativeSocketDispatcher
use, intrinsic :: iso_c_binding
implicit none


type, public :: FSICommNativeSocketDispatcher
     integer(kind=c_long_long )::reference
     contains
     procedure,public::createClientDispatcherInstanceForC
     procedure,public::createClientDispatcherInstance
     procedure,public::destroyDispatcherInstance
     
     
     	procedure,public::transferCoordinates
	procedure,private::transferCoordinates_internal
	procedure,public::transferData
	procedure,private::transferData_internal

end type FSICommNativeSocketDispatcher
contains
subroutine createClientDispatcherInstanceForC(this,host,port,buffer_size)
    class(FSICommNativeSocketDispatcher)::this
    character(kind=c_char),dimension(*)::host
    integer(kind=c_int)::port
    integer(kind=c_int)::buffer_size
    this%reference=0
    call fsi_fsicomm_f2c_nsd_create_client_instance(this%reference,host,port,buffer_size)
    
    

end subroutine createClientDispatcherInstanceForC

subroutine createClientDispatcherInstance(this,host,port,buffer_size)
    class(FSICommNativeSocketDispatcher)::this
    character(*)::host
    integer::port
    integer::buffer_size
    call this%createClientDispatcherInstanceForC(host//c_null_char,port,buffer_size)
    
    

end subroutine createClientDispatcherInstance

subroutine destroyDispatcherInstance(this)
     class(FSICommNativeSocketDispatcher)::this
     call fsi_fsicomm_f2c_nsd_destroy_instance(this%reference)

end subroutine destroyDispatcherInstance

subroutine transferData_internal(this,&
	data,data_len)
     use, intrinsic :: iso_c_binding
     class(FSICommNativeSocketDispatcher)::this
     	real(kind=c_double),intent(in),dimension(*)::data
	integer(kind=c_int),intent(in)::data_len

     call fsi_fsicomm_f2c_nsd_transferData(this%reference,&
data,data_len)
end subroutine transferData_internal

subroutine transferData(this,&
	data,data_len)
     use, intrinsic :: iso_c_binding
     class(FSICommNativeSocketDispatcher)::this
     	real(8),intent(in),dimension(*)::data
	integer,intent(in)::data_len

     call this%transferData_internal(&
data,data_len)
end subroutine transferData
subroutine transferCoordinates_internal(this,&
	coordId,coordId_len,&
	offsets,offsets_len,&
	hosts,hosts_len)
     use, intrinsic :: iso_c_binding
     class(FSICommNativeSocketDispatcher)::this
     	integer(kind=c_int),intent(in),dimension(*)::coordId
	integer(kind=c_int),intent(in)::coordId_len
	integer(kind=c_int),intent(in),dimension(*)::offsets
	integer(kind=c_int),intent(in)::offsets_len
	integer(kind=c_int),intent(in)::hosts_len
	type(c_ptr),dimension(*),intent(in)::hosts

     call fsi_fsicomm_f2c_nsd_transferCoordinates(this%reference,&
coordId,coordId_len,&
offsets,offsets_len,&
hosts,hosts_len)
end subroutine transferCoordinates_internal

subroutine transferCoordinates(this,&
	coordId,coordId_len,&
	offsets,offsets_len,&
	hosts,hosts_len)
     use, intrinsic :: iso_c_binding
     class(FSICommNativeSocketDispatcher)::this
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

     call this%transferCoordinates_internal(&
coordId,coordId_len,&
offsets,offsets_len,&
hostsPtrArray,hosts_len)
end subroutine transferCoordinates

end module  fsi_FSICommFNativeSocketDispatcher
