function(copy_target_property from to property)
    get_target_property(value ${from} ${property})
    set_target_properties(${to} PROPERTIES ${property} "${value}")
endfunction(copy_target_property)

