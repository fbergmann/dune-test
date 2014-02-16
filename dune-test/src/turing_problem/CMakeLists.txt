
set (HEADERS 

local_operator.hh   turing_initial.hh

)
set (SOURCES 
turing.cc
)

add_executable("dune_test" ${SOURCES})
target_link_dune_default_libraries("dune_test")

