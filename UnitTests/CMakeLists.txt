# File copied from group C
file(GLOB_RECURSE TEST_SRC "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/*.h")

add_executable(MolSimTests ${TEST_SRC})

target_link_libraries(
        MolSimTests
        PUBLIC
        GTest::gtest
        GTest::gtest_main
        MolSimLib # needs to be after gtest
)

target_include_directories(MolSimTests
        PRIVATE
        ${CMAKE_SOURCE_DIR}/src
)

include(GoogleTest)
gtest_discover_tests(MolSimTests)
