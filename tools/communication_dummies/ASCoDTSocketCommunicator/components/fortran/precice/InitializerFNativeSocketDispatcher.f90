module precice_InitializerFNativeSocketDispatcher
use, intrinsic :: iso_c_binding
implicit none


type, public :: InitializerNativeSocketDispatcher
     integer(kind=c_long_long )::reference
     contains
     procedure,public::createClientDispatcherInstanceForC
     procedure,public::createClientDispatcherInstance
     procedure,public::destroyDispatcherInstance
     
     
     	procedure,public::initializeAddresses
	procedure,private::initializeAddresses_internal
	procedure,public::initializeVertexes
	procedure,private::initializeVertexes_internal

end type InitializerNativeSocketDispatcher
contains
subroutine createClientDispatcherInstanceForC(this,host,port,buffer_size)
    class(InitializerNativeSocketDispatcher)::this
    character(kind=c_char),dimension(*)::host
    integer(kind=c_int)::port
    integer(kind=c_int)::buffer_size
    this%reference=0
    call precice_initializer_f2c_nsd_create_client_instance(this%reference,host,port,buffer_size)
    
    

end subroutine createClientDispatcherInstanceForC

subroutine createClientDispatcherInstance(this,host,port,buffer_size)
    class(InitializerNativeSocketDispatcher)::this
    character(*)::host
    integer::port
    integer::buffer_size
    call this%createClientDispatcherInstanceForC(host//c_null_char,port,buffer_size)
    
    

end subroutine createClientDispatcherInstance

subroutine destroyDispatcherInstance(this)
     class(InitializerNativeSocketDispatcher)::this
     call precice_initializer_f2c_nsd_destroy_instance(this%reference)

end subroutine destroyDispatcherInstance

subroutine initializeVertexes_internal(this,&
	vertexes,vertexes_len)
     use, intrinsic :: iso_c_binding
     class(InitializerNativeSocketDispatcher)::this
     	integer(kind=c_int),intent(in),dimension(*)::vertexes
	integer(kind=c_int),intent(in)::vertexes_len

     call precice_initializer_f2c_nsd_initializeVertexes(this%reference,&
vertexes,vertexes_len)
end subroutine initializeVertexes_internal

subroutine initializeVertexes(this,&
	vertexes,vertexes_len)
     use, intrinsic :: iso_c_binding
     class(InitializerNativeSocketDispatcher)::this
     	integer,intent(in),dimension(*)::vertexes
	integer,intent(in)::vertexes_len

     call this%initializeVertexes_internal(&
vertexes,vertexes_len)
end subroutine initializeVertexes
subroutine initializeAddresses_internal(this,&
	addresses,addresses_len)
     use, intrinsic :: iso_c_binding
     class(InitializerNativeSocketDispatcher)::this
     	integer(kind=c_int),intent(in)::addresses_len
	type(c_ptr),dimension(*),intent(in)::addresses

     call precice_initializer_f2c_nsd_initializeAddresses(this%reference,&
addresses,addresses_len)
end subroutine initializeAddresses_internal

subroutine initializeAddresses(this,&
	addresses,addresses_len)
     use, intrinsic :: iso_c_binding
     class(InitializerNativeSocketDispatcher)::this
     	character(*),intent(in),dimension(*)::addresses
	integer,intent(in)::addresses_len
	type(c_ptr),dimension(addresses_len) :: addressesPtrArray
	integer::addresses_ns
	character(255), dimension(addresses_len), target :: addressesFSArray
	do addresses_ns = 1, addresses_len
		addressesFSArray(addresses_ns) = addresses(addresses_ns)// C_NULL_CHAR
		addressesPtrArray(addresses_ns) = C_LOC(addressesFSArray(addresses_ns))
	end do

     call this%initializeAddresses_internal(&
addressesPtrArray,addresses_len)
end subroutine initializeAddresses

end module  precice_InitializerFNativeSocketDispatcher
