
set (HEADERS 

ccfv_local_operator.hh componentparameters.hh

)
set (SOURCES 
fish_ccfv.cc
)

add_executable(fish_ccfv ${SOURCES})
target_link_dune_default_libraries("fish_ccfv")

# copy config file if we don't have it
if (NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/fish_ccfv.conf")
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fish_ccfv.conf ${CMAKE_CURRENT_BINARY_DIR}/fish_ccfv.conf COPYONLY)
endif()

# create vtk dir, as otherwise the executable just falls over
if (NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/vtk")
file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/vtk")
endif()
