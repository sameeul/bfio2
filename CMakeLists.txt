cmake_minimum_required(VERSION 3.18)
project(libbfio2 LANGUAGES CXX)
#set(CMAKE_DEBUG_POSTFIX d)
set (CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION OFF)

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

set(CMAKE_CXX_FLAGS_RELEASE "-O3")

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake-modules")

#== Install OMEZarr PreReqs
option(WITH_OMEZARR "Build with OMEZarr Support" ON)
if(WITH_OMEZARR)
	set(ENV{Z5_INSTALL_DIR} "local_install")
	set(CMAKE_PREFIX_PATH ${CMAKE_BINARY_DIR}/$ENV{Z5_INSTALL_DIR})

	find_package(Boost)
	if(Boost_FOUND)
		set(ENV{Boost_INSTALLED} "TRUE")
	else()
		set(ENV{Boost_INSTALLED} "FALSE")
		message(STATUS "Boost will be downloaded and installed")
	endif()

	find_package(BLOSC)
	if(BLOSC_FOUND)
		set(ENV{BLOSC_INSTALLED} "TRUE")
	else()
		set(ENV{BLOSC_INSTALLED} "FALSE")
	endif()

	find_package(ZLIB)
	if(ZLIB_FOUND)
		set(ENV{ZLIB_INSTALLED} "TRUE")
	else()
		set(ENV{ZLIB_INSTALLED} "FALSE")
	endif()

	find_package(TIFF)
	if(TIFF_FOUND)
		set(ENV{TIFF_INSTALLED} "TRUE")
	else()
		set(ENV{TIFF_INSTALLED} "FALSE")
	endif()

	find_file(Z5_Installed "z5/z5.hxx")
	if(${Z5_Installed} STREQUAL "Z5_Installed-NOTFOUND")
		set(ENV{Z5_INSTALLED} "FALSE")
	else()
		set(ENV{Z5_INSTALLED} "TRUE")
		message(STATUS "Found Z5: ${Z5_Installed}")
	endif()

	if(UNIX)
		execute_process(COMMAND bash ${CMAKE_SOURCE_DIR}/build_scripts/install_prereq_linux.sh)
	endif()
	if(WIN32)
		execute_process(COMMAND cmd /C ${CMAKE_SOURCE_DIR}/build_scripts/install_prereq_win.bat)
	endif()

	find_package(Boost REQUIRED)
	if(Boost_FOUND)
		include_directories(${Boost_INCLUDE_DIR})
	endif()
	find_package(BLOSC REQUIRED)
	include_directories(${BLOSC_INCLUDE_DIR})
	list(APPEND bfio2_LIBRARIES ${BLOSC_LIBRARIES})
endif()

find_package(FastLoader REQUIRED)
if(FastLoader_FOUND)
	list(APPEND bfio2_LIBRARIES ${FastLoader_LIBRARIES})
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${FastLoader_CXX_FLAGS}")
	include_directories(${FastLoader_INCLUDE_DIR})
else()
	message(FATAL_ERROR "Unable to find FastLoader which is required for Image Loader.")
endif()
message(STATUS "CXX Flags" ${CMAKE_CXX_FLAGS})

#pugixml
include_directories(lib/pugixml)

#folly
include_directories(lib/folly)

#ThreadPool
include_directories(lib/ThreadPool)

set(SOURCE
	src/bfio2/ome_tiff_loader.h
	src/bfio2/ome_zarr_loader.h
	lib/pugixml/pugixml.cpp
	src/bfio2/utilities.cpp
)

# Python bindings.
if(BUILD_LIB)
	pybind11_add_module(libbfio2
		${SOURCE}
		src/bfio2/python/bfio2_interface.cpp
	)
endif()

# CLI
if(BUILD_CLI)
	add_executable(bfio2 ${SOURCE} src/bfio2/main.cpp)
	if(WITH_OMEZARR)
		target_compile_definitions(bfio2 PRIVATE -DWITH_BLOSC)
	endif()
	target_link_libraries(bfio2 PRIVATE stdc++fs)
	target_link_libraries(bfio2 PRIVATE ${bfio2_LIBRARIES})
	if (OPENMP_FOUND)
		target_link_libraries(bfio2 PRIVATE OpenMP::OpenMP_CXX)
	endif()
	if (APPLE)
		set_target_properties(bfio2 PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
	endif()
endif()

#Library
if(BUILD_LIB)
	if(WITH_OMEZARR)
		target_compile_definitions(libbfio2 PRIVATE -DWITH_BLOSC)
	endif()
	if (APPLE)
		set_target_properties(libbfio2 PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
	endif()
	target_compile_definitions(libbfio2 PRIVATE WITH_PYTHON_H)
	target_link_libraries(libbfio2 PRIVATE stdc++fs)
	target_link_libraries(libbfio2 PRIVATE ${bfio2_LIBRARIES})
	if (OPENMP_FOUND)
		target_link_libraries(libbfio2 PRIVATE OpenMP::OpenMP_CXX)
	endif()
endif()
