
set (HEADERS 

local_operator.hh componentparameters.hh initial_conditions.hh reactionadapter.hh

)
set (SOURCES 
turing_circle.cc
)

add_executable(turing_circle ${SOURCES})
target_link_dune_default_libraries("turing_circle")

# copy config file if we don't have it
if (NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/turing_circle.conf")
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/turing_circle.conf ${CMAKE_CURRENT_BINARY_DIR}/turing_circle.conf COPYONLY)
endif()

# create vtk dir, as otherwise the executable just falls over
if (NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/vtk")
file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/vtk")
endif()
