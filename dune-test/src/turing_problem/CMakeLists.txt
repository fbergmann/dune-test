
set (HEADERS 

local_operator.hh   turing_initial.hh

)
set (SOURCES 
turing.cc
)

add_executable("turing_problem" ${SOURCES})
target_link_dune_default_libraries("turing_problem")

# copy config file if we don't have it
if (NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/turing_problem.conf")
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/turing.conf ${CMAKE_CURRENT_BINARY_DIR}/turing_problem.conf COPYONLY)
endif()

# create vtk dir, as otherwise the executable just falls over
if (NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/vtk")
file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/vtk")
endif()
