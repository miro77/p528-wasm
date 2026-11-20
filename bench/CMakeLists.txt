###########################################
## BUILD THE PROFILE BENCHMARK
###########################################
set(BENCH_NAME "${LIB_NAME}Bench")
proplib_message("Configuring benchmark ${BENCH_NAME}")

add_executable(${BENCH_NAME} "P528Bench.cpp")
target_link_libraries(${BENCH_NAME} PRIVATE ${LIB_NAME})
configure_proplib_target(${BENCH_NAME})
if (MSVC)
    # the benchmark uses plain C file I/O for its reference files
    target_compile_definitions(${BENCH_NAME} PRIVATE _CRT_SECURE_NO_WARNINGS)
endif ()

# Run a small profile as part of the test suite. It fails if the profile API
# deviates from the point-by-point API, or if either of them deviates from the
# results stored before the profile optimizations were made.
if (RUN_TESTS)
    add_test(NAME ${BENCH_NAME}.ProfileMatchesReference
             COMMAND ${BENCH_NAME} --points 200 --check-reference
                     "${CMAKE_CURRENT_SOURCE_DIR}/reference_200.txt")
endif ()

proplib_message("Done configuring benchmark ${BENCH_NAME}")
