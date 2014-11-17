module precice_ReceiverFNativeSocketDispatcher
use, intrinsic :: iso_c_binding
implicit none


type, public :: ReceiverNativeSocketDispatcher
     integer(kind=c_long_long )::reference
     contains
     procedure,public::createClientDispatcherInstanceForC
     procedure,public::createClientDispatcherInstance
     procedure,public::destroyDispatcherInstance
     
     
     	procedure,public::receive
	procedure,private::receive_internal

end type ReceiverNativeSocketDispatcher
contains
subroutine createClientDispatcherInstanceForC(this,host,port,buffer_size)
    class(ReceiverNativeSocketDispatcher)::this
    character(kind=c_char),dimension(*)::host
    integer(kind=c_int)::port
    integer(kind=c_int)::buffer_size
    this%reference=0
    call precice_receiver_f2c_nsd_create_client_instance(this%reference,host,port,buffer_size)
    
    

end subroutine createClientDispatcherInstanceForC

subroutine createClientDispatcherInstance(this,host,port,buffer_size)
    class(ReceiverNativeSocketDispatcher)::this
    character(*)::host
    integer::port
    integer::buffer_size
    call this%createClientDispatcherInstanceForC(host//c_null_char,port,buffer_size)
    
    

end subroutine createClientDispatcherInstance

subroutine destroyDispatcherInstance(this)
     class(ReceiverNativeSocketDispatcher)::this
     call precice_receiver_f2c_nsd_destroy_instance(this%reference)

end subroutine destroyDispatcherInstance

subroutine receive_internal(this,&
	data,&
	index,&
	rank,&
	tag)
     use, intrinsic :: iso_c_binding
     class(ReceiverNativeSocketDispatcher)::this
     	real(kind=c_double),intent(in)::data
	integer(kind=c_int),intent(in)::index
	integer(kind=c_int),intent(in)::rank
	integer(kind=c_int),intent(inout)::tag

     call precice_receiver_f2c_nsd_receive(this%reference,&
data,&
index,&
rank,&
tag)
end subroutine receive_internal

subroutine receive(this,&
	data,&
	index,&
	rank,&
	tag)
     use, intrinsic :: iso_c_binding
     class(ReceiverNativeSocketDispatcher)::this
     	real(8),intent(in)::data
	integer,intent(in)::index
	integer,intent(in)::rank
	integer,intent(inout)::tag

     call this%receive_internal(&
data,&
index,&
rank,&
tag)
end subroutine receive

end module  precice_ReceiverFNativeSocketDispatcher
