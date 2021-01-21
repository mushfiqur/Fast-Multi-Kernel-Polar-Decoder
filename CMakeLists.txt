# CMakeList.txt : CMake project for SoftwareModel, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

project (SoftwareModel)

# Add source to this project's executable.

file(GLOB SOURCES "src/*.cpp")
add_executable (SoftwareModel ${SOURCES})

# set(CMAKE_CXX_FLAGS "-O3")

target_include_directories(SoftwareModel
	PRIVATE
		${CMAKE_CURRENT_SOURCE_DIR}/inc
		${CMAKE_CURRENT_SOURCE_DIR}/inc/Decoders
		${CMAKE_CURRENT_SOURCE_DIR}/inc/Utils
		${CMAKE_CURRENT_SOURCE_DIR}/inc/Utils/Algo
		${CMAKE_CURRENT_SOURCE_DIR}/inc/Eigen
)

# TODO: Add tests and install targets if needed.
