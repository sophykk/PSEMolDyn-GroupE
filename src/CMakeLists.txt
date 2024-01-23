# collect all cpp files
file(GLOB_RECURSE MY_SRC
        "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/*.cxx"
        # header don't need to be included but this might be necessary for some IDEs
        "${CMAKE_CURRENT_SOURCE_DIR}/*.h"
        EXCLUDE  "${CMAKE_CURRENT_SOURCE_DIR}/MolSim.cpp"
)

# Make a MolSim library so we don't need to recompile everything when testing
add_library(MolSimLib ${MY_SRC})

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

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
        OpenMP::OpenMP_CXX
)

# create make target
add_executable(MolSim MolSim.cpp)

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


# enable OpenMP if available
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(MolSim PUBLIC OpenMP::OpenMP_CXX)
endif()

# activate all compiler warnings. Clean up your code :P
# depending on the compiler different flags are used
target_compile_options(MolSim
        PRIVATE
        $<$<CXX_COMPILER_ID:GNU>:-Wall>
        # disable some spam warnings for icpc...
        $<$<CXX_COMPILER_ID:Intel>:-w3 -wd383,981,1418,1572,2259>
)
