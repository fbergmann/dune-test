
set (HEADERS 

local_operator.hh componentparameters.hh initial_conditions.hh reactionadapter.hh

)
set (SOURCES 
pnas_fish.cc
)

add_executable(pnas_fish ${SOURCES})
target_link_dune_default_libraries("pnas_fish")

# copy config file if we don't have it
if (NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/pnas_fish.conf")
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/pnas_fish.conf ${CMAKE_CURRENT_BINARY_DIR}/pnas_fish.conf COPYONLY)
endif()

# create vtk dir, as otherwise the executable just falls over
if (NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/vtk")
file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/vtk")
endif()
