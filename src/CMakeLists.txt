# Collect every .cxx file in src/ and its sub-dirs
file(GLOB_RECURSE FRAMEWORK_SRC CONFIGURE_DEPENDS
     "${CMAKE_CURRENT_SOURCE_DIR}/*.cxx")

# Keep headers visible to clients
set(FRAMEWORK_INC_DIR ${PROJECT_SOURCE_DIR}/include)
file(GLOB_RECURSE FRAMEWORK_HEADERS CONFIGURE_DEPENDS
     "${FRAMEWORK_INC_DIR}/*.hxx" "${FRAMEWORK_INC_DIR}/*.hpp")

add_library(AnalysisModules SHARED ${FRAMEWORK_SRC} ${FRAMEWORK_HEADERS})

target_include_directories(AnalysisModules
    PUBLIC
      ${FRAMEWORK_INC_DIR}
)

target_link_libraries(AnalysisModules
    PUBLIC
      ROOT::Core ROOT::RIO ROOT::Tree
      ROOT::Hist ROOT::Imt ROOT::ROOTDataFrame
      ROOT::TMVA ROOT::RooFit ROOT::RooStats
)

# Optional dictionary if you ever need I/O for a custom class ----------------
# set(DICT_HEADERS ${FRAMEWORK_INC_DIR}/Framework/MyClass.hxx)
# root_generate_dictionary(G__bsm_dict  "${DICT_HEADERS}"
#                          LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/LinkDef.h)
# target_sources(AnalysisModules PRIVATE ${G__bsm_dict})
# target_link_libraries(AnalysisModules PRIVATE ROOT::Core)