cmake_minimum_required(VERSION 2.8)
project(SibeliaZ)
set(CMAKE_BUILD_TYPE RELEASE)
set(CMAKE_CXX_STANDARD 11)

if (NOT NO_DEPENDENCIES)
	option(spoa_build_executable "Build spoa standalone tool" ON)
	set(spoa_build_executable ON)
	add_subdirectory(spoa)
	add_subdirectory(TwoPaCo/src)
endif()

add_subdirectory(SibeliaZ-LCB)
include(CPack)
