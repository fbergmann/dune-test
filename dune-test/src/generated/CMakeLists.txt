
set (HEADERS 

local_operator.hh componentparameters.hh initial_conditions.hh reactionadapter.hh

)
set (SOURCES 
main.cc
)

add_executable(main ${SOURCES})
target_link_dune_default_libraries("main")

# copy config file if we don't have it
if (NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/main.conf")
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/main.conf ${CMAKE_CURRENT_BINARY_DIR}/main.conf COPYONLY)
endif()

# create vtk dir, as otherwise the executable just falls over
if (NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/vtk")
file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/vtk")
endif()
