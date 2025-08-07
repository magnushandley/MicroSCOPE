# helper that turns a 1-line .cxx into an executable
function(make_runner target)
  add_executable(${target} ${target}.cxx)

  # <── THIS WAS MISSING ────────────────────────────────────────────────┐
  target_include_directories(${target} PRIVATE
      ${PROJECT_SOURCE_DIR}/include   # Framework and module headers
      ${ROOT_INCLUDE_DIRS}            # (optional but nice)
  )                                   #┘

  target_link_libraries(${target} PRIVATE AnalysisModules)
endfunction()

make_runner(run_slimmer)
make_runner(run_preselection)

#add_executable(run_slimmer run_slimmer.cxx)

#target_include_directories(run_slimmer PRIVATE
#    ${PROJECT_SOURCE_DIR}/include   # Framework and module headers
#    ${ROOT_INCLUDE_DIRS}            # (optional but nice)
#)

#target_link_libraries(run_slimmer PRIVATE AnalysisModules)
