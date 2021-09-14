add_executable(parser_test parser_test.cc)
target_link_libraries(parser_test PUBLIC gtest_main parser)

add_test(
    NAME parser_test
    COMMAND parser_test
)

add_executable(pulse_test pulse_test.cc)
target_link_libraries(pulse_test PUBLIC gtest_main pulse)

add_test(
    NAME pulse_test
    COMMAND pulse_test
)
