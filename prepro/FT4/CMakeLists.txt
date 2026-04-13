file(GLOB FT4_SRCS "*.cpp")
add_executable(ft4 ft4.cpp)
target_include_directories(ft4 PUBLIC ${COMMON_HEADERS})
target_link_libraries(ft4 librock)
install(TARGETS ft4 DESTINATION bin)

