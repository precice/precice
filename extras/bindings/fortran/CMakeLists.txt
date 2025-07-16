set_property(TARGET precice PROPERTY PUBLIC_HEADER
  ${CMAKE_CURRENT_LIST_DIR}/include/precice/preciceFortran.hpp
  APPEND)

target_sources(preciceCore PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/src/preciceFortran.cpp
  ${CMAKE_CURRENT_LIST_DIR}/include/precice/preciceFortran.hpp
  )

target_include_directories(preciceCore PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>
  )

target_include_directories(precice PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>
  )


# Mock

set_property(TARGET preciceMocked PROPERTY PUBLIC_HEADER
  ${CMAKE_CURRENT_LIST_DIR}/include/precice/preciceFortran.hpp
  APPEND)

target_sources(preciceMocked PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/include/precice/preciceFortran.hpp
  ${CMAKE_CURRENT_LIST_DIR}/src/preciceFortran.cpp
  )

target_include_directories(preciceMocked PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>
  )
