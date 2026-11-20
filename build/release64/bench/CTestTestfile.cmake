# CMake generated Testfile for 
# Source directory: M:/test_repos/p528/bench
# Build directory: M:/test_repos/p528/build/release64/bench
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
if(CTEST_CONFIGURATION_TYPE MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
  add_test(P528Bench.ProfileMatchesReference "M:/test_repos/p528/bin/Release/P528Bench.exe" "--points" "200" "--check-reference" "M:/test_repos/p528/bench/reference_200.txt")
  set_tests_properties(P528Bench.ProfileMatchesReference PROPERTIES  _BACKTRACE_TRIPLES "M:/test_repos/p528/bench/CMakeLists.txt;19;add_test;M:/test_repos/p528/bench/CMakeLists.txt;0;")
else()
  add_test(P528Bench.ProfileMatchesReference NOT_AVAILABLE)
endif()
