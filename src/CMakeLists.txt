# collect all cpp files
file(GLOB_RECURSE MY_SRC
        "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
        # header don't need to be included but this might be necessary for some IDEs
        "${CMAKE_CURRENT_SOURCE_DIR}/*.h"
        EXCLUDE  "${CMAKE_CURRENT_SOURCE_DIR}/MolSim.cpp"
)

# Make a MolSim library so we don't need to recompile everything when testing
add_library(MolSimLib ${MY_SRC})

target_include_directories(MolSimLib
        PUBLIC
        ${XercesC_INCLUDE_DIRS}
        ${CMAKE_SOURCE_DIR}/libs/libxsd
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(MolSimLib
        # stuff that is used in headers and source files
        PUBLIC
        ${XercesC_LIBRARY}
        spdlog::spdlog
)

# create make target
add_executable(MolSim MolSim.cpp)

# set cxx standard. You may raise this if you want.
target_compile_features(MolSim
        PRIVATE
        cxx_std_17
)

target_include_directories(MolSim
        PUBLIC
        ${CMAKE_SOURCE_DIR}/libs/libxsd
        ${spdlog_SOURCE_DIR}/include/spdlog
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(MolSim
        PUBLIC
        MolSimLib
)

# activate all compiler warnings. Clean up your code :P
# depending on the compiler different flags are used
target_compile_options(MolSim
        PRIVATE
        $<$<CXX_COMPILER_ID:GNU>:-Wall>
        # disable some spam warnings for icpc...
        $<$<CXX_COMPILER_ID:Intel>:-w3 -wd383,981,1418,1572,2259>
)
