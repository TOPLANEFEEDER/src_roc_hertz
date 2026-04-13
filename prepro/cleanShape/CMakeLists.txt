file(GLOB cleanShape_SRCS "*.cpp")
add_executable(cleanShape cleanShape.cpp)
target_include_directories(cleanShape PUBLIC ${COMMON_HEADERS})
target_link_libraries(cleanShape librock)
install(TARGETS cleanShape DESTINATION bin)

