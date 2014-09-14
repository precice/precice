module precice_CommunicatorFNativeSocketDispatcher
use, intrinsic :: iso_c_binding
implicit none


type, public :: CommunicatorNativeSocketDispatcher
     integer(kind=c_long_long )::reference
     contains
     procedure,public::createClientDispatcherInstanceForC
     procedure,public::createClientDispatcherInstance
     procedure,public::destroyDispatcherInstance
     
     
     	procedure,public::setData
	procedure,private::setData_internal

end type CommunicatorNativeSocketDispatcher
contains
subroutine createClientDispatcherInstanceForC(this,host,port,buffer_size)
    class(CommunicatorNativeSocketDispatcher)::this
    character(kind=c_char),dimension(*)::host
    integer(kind=c_int)::port
    integer(kind=c_int)::buffer_size
    this%reference=0
    call precice_communicator_f2c_nsd_create_client_instance(this%reference,host,port,buffer_size)
    
    

end subroutine createClientDispatcherInstanceForC

subroutine createClientDispatcherInstance(this,host,port,buffer_size)
    class(CommunicatorNativeSocketDispatcher)::this
    character(*)::host
    integer::port
    integer::buffer_size
    call this%createClientDispatcherInstanceForC(host//c_null_char,port,buffer_size)
    
    

end subroutine createClientDispatcherInstance

subroutine destroyDispatcherInstance(this)
     class(CommunicatorNativeSocketDispatcher)::this
     call precice_communicator_f2c_nsd_destroy_instance(this%reference)

end subroutine destroyDispatcherInstance

subroutine setData_internal(this,&
	data,data_len)
     use, intrinsic :: iso_c_binding
     class(CommunicatorNativeSocketDispatcher)::this
     	real(kind=c_double),intent(in),dimension(*)::data
	integer(kind=c_int),intent(in)::data_len

     call precice_communicator_f2c_nsd_setData(this%reference,&
data,data_len)
end subroutine setData_internal

subroutine setData(this,&
	data,data_len)
     use, intrinsic :: iso_c_binding
     class(CommunicatorNativeSocketDispatcher)::this
     	real(8),intent(in),dimension(*)::data
	integer,intent(in)::data_len

     call this%setData_internal(&
data,data_len)
end subroutine setData

end module  precice_CommunicatorFNativeSocketDispatcher
