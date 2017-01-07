
cmake_minimum_required (VERSION 2.6)
project (oshun_c_routines_for_speedup)

#FIND_PACKAGE(CUDA)
find_program(PYTHON "python")

# if this is windows, then define NUMPY1_EXPORTS to make the shared library...
if( MSVC )
	add_definitions("/DNUMPY1_EXPORTS")
endif( MSVC )

# Set the default build type to 'Release' (fastest) and amke a nice string explaining hte different build-types to the user
# If they are using a gui to configure.
if(NOT CMAKE_BUILD_TYPE)
       set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE)

# Zip up the python code...
file(TO_NATIVE_PATH "${CMAKE_SOURCE_DIR}/utils/build/ZipUpModule.py" ZIP_BUILDER)
file(TO_NATIVE_PATH "${CMAKE_SOURCE_DIR}/source/python" PYTHON_SOURCE)
file(TO_NATIVE_PATH "${CMAKE_SOURCE_DIR}/source/python/oshun_modules" OSHUN_PYTHON_SOURCE)
file(TO_NATIVE_PATH "${CMAKE_SOURCE_DIR}/utils/build/AddFilesToZip.py" ADD_TO_ZIP)


add_custom_target(	oshun ALL 
					"${PYTHON}" "${ZIP_BUILDER}" "${PYTHON_SOURCE}" "oshun" "oshun_modules" "${OSHUN_PYTHON_SOURCE}"
					COMMENT "Zipping contents of /source/python/oshun_modules into package oshun.zip" )

add_library( oshun_c_routines_for_speedup SHARED "${CMAKE_SOURCE_DIR}/source/c/numpy1.cpp")

set_target_properties( oshun_c_routines_for_speedup PROPERTIES OUTPUT_NAME oshun_c_routines_for_speedup PREFIX "")

add_custom_target(	add_c_to_oshun ALL
					"${PYTHON}" "${ADD_TO_ZIP}" "oshun.zip" $<TARGET_FILE:oshun_c_routines_for_speedup> "oshun_modules/shared_libraries"
					COMMENT "Injecting the Shared Library Into the oshun.zip package" )

ADD_DEPENDENCIES(oshun_c_routines_for_speedup  oshun )
ADD_DEPENDENCIES(add_c_to_oshun oshun_c_routines_for_speedup  )

#MESSAGE( "Cmake reports a compilier with ID ${CMAKE_CXX_COMPILER_ID}" )

if (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
	set_source_files_properties(numpy1.cpp PROPERTIES COMPILE_FLAGS -fno-alias)
	#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}  -fno-alias")

elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
	set_source_files_properties(numpy1.cpp PROPERTIES COMPILE_FLAGS -fstrict-aliasing)
	#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}  -fstrict-aliasing")
endif()

#MESSAGE( "Release flags for c++: ${CMAKE_CXX_FLAGS_RELEASE}" )
#MESSAGE( "  Debug flags for c++: ${CMAKE_CXX_FLAGS_DEBUG}" )

#file(RENAME "${CMAKE_BINARY_DIR}/libnumpy1.so" "${CMAKE_BINARY_DIR}/numpy1.so")
