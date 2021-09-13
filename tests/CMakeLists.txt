add_executable(Ptrac_test MCNPGeometryPtrac_test.cc)
target_link_libraries(Ptrac_test PUBLIC gtest_main MCNPGeometry)

add_executable(parser_test parser_test.cc)
target_link_libraries(parser_test PUBLIC gtest_main parser)


add_test(
    NAME Ptrac_test
    COMMAND Ptrac_test
)

add_test(
    NAME parser_test
    COMMAND parser_test
)