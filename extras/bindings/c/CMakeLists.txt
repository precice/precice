set_property(TARGET precice PROPERTY PUBLIC_HEADER
  ${CMAKE_CURRENT_LIST_DIR}/include/precice/preciceC.h
  APPEND)

target_sources(preciceCore PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/include/precice/preciceC.h
  ${CMAKE_CURRENT_LIST_DIR}/src/preciceC.cpp
  )

target_include_directories(preciceCore PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>
  )

target_include_directories(precice PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>
  )
