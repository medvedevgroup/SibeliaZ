project(graphdump CXX)
cmake_minimum_required(VERSION 2.8)
set(CMAKE_BUILD_TYPE RELEASE)
set(CMAKE_PROJECT_NAME sibeliaz)
set(CMAKE_CXX_STANDARD 11)


set(twopaco_SOURCE_DIR ../TwoPaCo/src/common)
add_executable(sibeliaz-lcb sibeliaz.cpp blocksfinder.cpp ${twopaco_SOURCE_DIR}/dnachar.cpp ${twopaco_SOURCE_DIR}/streamfastaparser.cpp)
include_directories(${twopaco_SOURCE_DIR})
find_package(OpenMP)
target_link_libraries(sibeliaz-lcb PUBLIC OpenMP::OpenMP_CXX)
install(TARGETS sibeliaz-lcb RUNTIME DESTINATION bin)
install(PROGRAMS sibeliaz DESTINATION bin)
