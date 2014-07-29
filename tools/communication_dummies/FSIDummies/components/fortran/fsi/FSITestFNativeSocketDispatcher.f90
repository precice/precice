module fsi_FSITestFNativeSocketDispatcher
use, intrinsic :: iso_c_binding
implicit none


type, public :: FSITestNativeSocketDispatcher
     integer(kind=c_long_long )::reference
     contains
     procedure,public::createClientDispatcherInstanceForC
     procedure,public::createClientDispatcherInstance
     procedure,public::destroyDispatcherInstance
     
     
     	procedure,public::test
	procedure,private::test_internal

end type FSITestNativeSocketDispatcher
contains
subroutine createClientDispatcherInstanceForC(this,host,port,buffer_size)
    class(FSITestNativeSocketDispatcher)::this
    character(kind=c_char),dimension(*)::host
    integer(kind=c_int)::port
    integer(kind=c_int)::buffer_size
    this%reference=0
    call fsi_fsitest_f2c_nsd_create_client_instance(this%reference,host,port,buffer_size)
    
    

end subroutine createClientDispatcherInstanceForC

subroutine createClientDispatcherInstance(this,host,port,buffer_size)
    class(FSITestNativeSocketDispatcher)::this
    character(*)::host
    integer::port
    integer::buffer_size
    call this%createClientDispatcherInstanceForC(host//c_null_char,port,buffer_size)
    
    

end subroutine createClientDispatcherInstance

subroutine destroyDispatcherInstance(this)
     class(FSITestNativeSocketDispatcher)::this
     call fsi_fsitest_f2c_nsd_destroy_instance(this%reference)

end subroutine destroyDispatcherInstance

subroutine test_internal(this)
     use, intrinsic :: iso_c_binding
     class(FSITestNativeSocketDispatcher)::this
     
     call fsi_fsitest_f2c_nsd_test(this%reference)
end subroutine test_internal

subroutine test(this)
     use, intrinsic :: iso_c_binding
     class(FSITestNativeSocketDispatcher)::this
     
     call this%test_internal()
end subroutine test

end module  fsi_FSITestFNativeSocketDispatcher
