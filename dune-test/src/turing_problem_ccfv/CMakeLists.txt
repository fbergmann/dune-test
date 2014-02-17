
set (HEADERS 

ccfv_local_operator.hh componentparameters.hh

)
set (SOURCES 
turing_ccfv.cc
)

add_executable(turing_ccfv ${SOURCES})
target_link_dune_default_libraries("turing_ccfv")

# copy config file if we don't have it
if (NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/turing_ccfv.conf")
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/turing_ccfv.conf ${CMAKE_CURRENT_BINARY_DIR}/turing_ccfv.conf COPYONLY)
endif()