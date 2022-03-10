cmake_minimum_required(VERSION 3.18)
project(libbfio2 LANGUAGES CXX)
set(CMAKE_DEBUG_POSTFIX d)
set (CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#OpenMP
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
	message(STATUS "OpenMP Found: ${OpenMP_CXX_FLAGS}")
else()
	message(STATUS "OpenMP Not Found")
endif()

# PyBind
add_subdirectory(lib/pybind11)
include_directories(
	lib/pybind11/include/
)

# Release Type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS_RELEASE "-O2")

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake-modules")

find_package(FastLoader REQUIRED)
if(FastLoader_FOUND)
	list(APPEND bfio2_LIBRARIES ${FastLoader_LIBRARIES})
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${FastLoader_CXX_FLAGS}")
	include_directories(${FastLoader_INCLUDE_DIR})
else()
	message(FATAL_ERROR "Unable to find FastLoader which is required for Image Loader.")
endif()
message(STATUS ${CMAKE_CXX_FLAGS})

#pugixml
include_directories(lib/pugixml)


set(SOURCE 
	src/bfio2/ome_tiff_loader.cpp 
	lib/pugixml/pugixml.cpp
)

# Python bindings.
pybind11_add_module(libbfio2 
	${SOURCE}
	src/bfio2/python/bfio2_interface.cpp
)

# CLI
if(BUILD_CLI)
	add_executable(bfio2 ${SOURCE} src/bfio2/main.cpp)
	target_link_libraries(bfio2 PRIVATE stdc++fs)
	target_link_libraries(bfio2 PRIVATE ${bfio2_LIBRARIES})	
endif()

target_link_libraries(libbfio2 PRIVATE stdc++fs)
target_link_libraries(libbfio2 PRIVATE ${bfio2_LIBRARIES})

