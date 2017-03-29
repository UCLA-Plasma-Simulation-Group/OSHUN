cmake_minimum_required(VERSION 2.8)
project(oshun1d)

file(COPY input/inputdeck DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

add_compile_options(-std=c++0x)

find_package(MPI REQUIRED COMPONENTS CXX)

# HDF5 needs to be compiled with C++ enabled
# E.g. - On hoffman:
# 'cd $OSHUN'
# 'mkdir $buildname'
# 'module load hdf5'
# 'export CXX=MPICXX'
# 'cmake ..'
# 'make'
#        source/laser.cpp
#source/laser.h

find_package(HDF5 REQUIRED COMPONENTS CXX)
include_directories(${HDF5_INCLUDE_DIRS})

set(SOURCE_FILES
        source/collisions.cpp
        source/collisions.h
        source/export.cpp
        source/export.h
        source/exprtk.hpp
        source/formulary.cpp
        source/formulary.h
        source/implicitE.cpp
        source/implicitE.h
        source/input.cpp
        source/input.h
#        source/interspeciescollisions.cpp
#        source/interspeciescollisions.h
        source/fluid.cpp
        source/fluid.h
        source/lib-algorithms.h
        source/lib-array.h
        source/main.cpp
        source/nmethods.cpp
        source/nmethods.h
        source/parallel.cpp
        source/parallel.h
        source/setup.cpp
        source/setup.h
        source/state.cpp
        source/state.h
        source/vlasov.cpp
        source/vlasov.h source/functors.cpp source/functors.h)

add_executable(oshun1d ${SOURCE_FILES})

## Compiler flags
if(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "-O2")        ## Optimize

endif()

target_link_libraries(oshun1d ${MPI_LIBRARIES} ${HDF5_CXX_LIBRARIES})
