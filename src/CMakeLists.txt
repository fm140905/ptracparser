add_library(parser SHARED parser.cc)

add_library(pulse SHARED pulse.cc)
target_link_libraries(pulse PUBLIC parser)

add_executable(main main.cc)
target_link_libraries(main PUBLIC parser pulse)