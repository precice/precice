
# This function overrides the option OPTION with the value of the TPL_OPTION if it is defined.
function(xsdk_tpl_option_override OPTION TPL_OPTION)
  if(DEFINED ${TPL_OPTION})
    set(${OPTION} ${${TPL_OPTION}} CACHE BOOL "TPL override!" FORCE)
  endif()
endfunction()

# This function requires a set of variables for a given TPL module of NAME.
# It will first list the state of every variable and then throw an error.
function(xsdk_tpl_require NAME)
  set(TPL_REQ_MISSING "")

  message(STATUS "Checking requirements for TPL override of ${NAME}")

  # Lists all required variables alongside their state.
  foreach(TPL_REQ ${ARGN})
    if ( NOT DEFINED ${TPL_REQ} )
      message(STATUS "  * ${TPL_REQ} - NOT DEFINED")
      list(APPEND TPL_REQ_MISSING "${TPL_REQ}")
    else()
      message(STATUS "  * ${TPL_REQ} - defined")
    endif()
  endforeach(TPL_REQ)

  if(TPL_REQ_MISSING)
    message(FATAL_ERROR "Checking requirements for TPL override of ${NAME} - failure")
  else()
    message(STATUS "Checking requirements for TPL override of ${NAME} - success")
  endif()
endfunction()
