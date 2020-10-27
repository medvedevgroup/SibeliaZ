cmake_minimum_required(VERSION 2.8)
set(CMAKE_POLICY_DEFAULT_CMP0048 NEW)

include(ExternalProject)


set(CMAKE_POLICY_WARNING_CMP0048 OFF)
set(CMAKE_BUILD_TYPE RELEASE)
set(CMAKE_CXX_STANDARD 11)

if (NOT NO_DEPENDENCIES)
	option(spoa_build_executable "Build spoa standalone tool" ON)
	set(spoa_build_executable ON)
	add_subdirectory(spoa)
	add_subdirectory(TwoPaCo/src)
	ExternalProject_Add(maf2synteny SOURCE_DIR ${PROJECT_SOURCE_DIR}/maf2synteny CONFIGURE_COMMAND "" BUILD_COMMAND make BUILD_IN_SOURCE 1 INSTALL_COMMAND "")
	install(PROGRAMS ${PROJECT_SOURCE_DIR}/maf2synteny/maf2synteny DESTINATION bin)
endif()

add_subdirectory(SibeliaZ-LCB)
include(CPack)
