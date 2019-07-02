# Copies a property from a target to another
function(copy_target_property from to property)
  get_target_property(value ${from} ${property})
  if(value)
    set_target_properties(${to} PROPERTIES ${property} "${value}")
  endif(value)
endfunction(copy_target_property)

