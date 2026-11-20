############################################
## CONFIGURE UNIT TESTS
############################################
set(TEST_NAME "${LIB_NAME}Test")
proplib_message("Configuring library tests ${TEST_NAME}")

add_executable(
    ${TEST_NAME}
    P528GTest.h
    P528GTestUtils.cpp
    TestP528.cpp
    TestP528DataTables.cpp
    TestP528Profile.cpp
)


# Set PropLib compiler option defaults
configure_proplib_target(${TEST_NAME})

###########################################
## SET UP AND DISCOVER TESTS
###########################################
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

# Explicitly link the filesystem library (only exists/needed for GCC's libstdc++)
target_link_libraries(${TEST_NAME} ${LIB_NAME} GTest::gtest_main
    "$<$<CXX_COMPILER_ID:GNU>:stdc++fs>")
include(GoogleTest)
gtest_discover_tests(${TEST_NAME})

proplib_message("Done configuring library tests ${TEST_NAME}")
