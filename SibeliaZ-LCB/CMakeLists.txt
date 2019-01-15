project(graphdump CXX)
cmake_minimum_required(VERSION 2.8)
set(CMAKE_BUILD_TYPE RELEASE)
set(CMAKE_PROJECT_NAME sibelia-pp)
if(CMAKE_COMPILER_IS_GNUCXX)
	list(APPEND "CMAKE_CXX_FLAGS" "-std=c++0x")
endif()

set(twopaco_SOURCE_DIR ../TwoPaCo/src/common)
add_executable(sibeliaz-lcb sibeliaz.cpp blocksfinder.cpp ${twopaco_SOURCE_DIR}/dnachar.cpp ${twopaco_SOURCE_DIR}/streamfastaparser.cpp)
link_directories(${TBB_LIB_DIR})
include_directories(${twopaco_SOURCE_DIR} ${TBB_INCLUDE_DIR})
target_link_libraries(sibeliaz-lcb "tbb")

