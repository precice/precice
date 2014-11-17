module precice_MainFNativeSocketDispatcher
use, intrinsic :: iso_c_binding
implicit none


type, public :: MainNativeSocketDispatcher
     integer(kind=c_long_long )::reference
     contains
     procedure,public::createClientDispatcherInstanceForC
     procedure,public::createClientDispatcherInstance
     procedure,public::destroyDispatcherInstance
     
     
     	procedure,public::main
	procedure,private::main_internal

end type MainNativeSocketDispatcher
contains
subroutine createClientDispatcherInstanceForC(this,host,port,buffer_size)
    class(MainNativeSocketDispatcher)::this
    character(kind=c_char),dimension(*)::host
    integer(kind=c_int)::port
    integer(kind=c_int)::buffer_size
    this%reference=0
    call precice_main_f2c_nsd_create_client_instance(this%reference,host,port,buffer_size)
    
    

end subroutine createClientDispatcherInstanceForC

subroutine createClientDispatcherInstance(this,host,port,buffer_size)
    class(MainNativeSocketDispatcher)::this
    character(*)::host
    integer::port
    integer::buffer_size
    call this%createClientDispatcherInstanceForC(host//c_null_char,port,buffer_size)
    
    

end subroutine createClientDispatcherInstance

subroutine destroyDispatcherInstance(this)
     class(MainNativeSocketDispatcher)::this
     call precice_main_f2c_nsd_destroy_instance(this%reference)

end subroutine destroyDispatcherInstance

subroutine main_internal(this)
     use, intrinsic :: iso_c_binding
     class(MainNativeSocketDispatcher)::this
     
     call precice_main_f2c_nsd_main(this%reference)
end subroutine main_internal

subroutine main(this)
     use, intrinsic :: iso_c_binding
     class(MainNativeSocketDispatcher)::this
     
     call this%main_internal()
end subroutine main

end module  precice_MainFNativeSocketDispatcher
