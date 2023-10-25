add_library(parser STATIC parser.cc)

add_library(pulse STATIC pulse.cc)
target_link_libraries(pulse PUBLIC parser)

add_executable(main main.cc)
target_link_libraries(main PUBLIC parser pulse)
# set_target_properties(main PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")