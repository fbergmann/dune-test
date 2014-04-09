add_subdirectory(turing_problem)
add_subdirectory(turing_problem_ccfv)
add_subdirectory(fish_ccfv)

if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/generated)
  add_subdirectory(generated)
endif()
