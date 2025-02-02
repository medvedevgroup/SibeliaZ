project(SibeliaZ-LCB CXX)
set(CMAKE_BUILD_TYPE RELEASE)
set(CMAKE_PROJECT_NAME sibeliaz)

set(common_SOURCE_DIR ./common)
add_executable(sibeliaz-lcb sibeliaz.cpp blocksfinder.cpp ${common_SOURCE_DIR}/dnachar.cpp ${common_SOURCE_DIR}/streamfastaparser.cpp)

include_directories(${common_SOURCE_DIR} ${TBB_LIB_DIR})
find_package(OpenMP)
target_link_libraries(sibeliaz-lcb PUBLIC OpenMP::OpenMP_CXX)
install(TARGETS sibeliaz-lcb RUNTIME DESTINATION bin)
install(PROGRAMS sibeliaz DESTINATION bin)
