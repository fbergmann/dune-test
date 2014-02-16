
set (HEADERS 

ccfv_local_operator.hh componentparameters

)
set (SOURCES 
turing.cc
)

add_executable("dune_test" ${SOURCES})
target_link_dune_default_libraries("dune_test")

