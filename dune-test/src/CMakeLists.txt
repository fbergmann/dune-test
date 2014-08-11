
foreach(dir 
        turing_circle
        turing_current
        turing_problem
        turing_problem_ccfv
        fish_ccfv
        generated2
        Jana_wolf
        GepasiYeast
        pnas_fish
        nfkb
        nfkb2
)

  if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${dir})
    add_subdirectory(${dir})
  endif()

endforeach()
