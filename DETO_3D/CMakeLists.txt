cmake_minimum_required(VERSION 3.0)

project(DETO VERSION 1.0 LANGUAGES CXX)

find_package(MPI REQUIRED)
if (MPI_FOUND)
  set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
endif()

set(SOURCES 
  src/main.cpp
  src/deto.cpp
  src/universe.cpp
  src/inputdeto.cpp
  src/output.cpp
  src/lammpsIO.cpp
  src/optimize.cpp
  src/simulations.cpp
  src/error.cpp)
  
find_library(LAMMPS_LIBRARIES lammps PATHS ${CMAKE_SOURCE_DIR}/lammps/build NO_DEFAULT_PATH)


add_executable(deto ${SOURCES})
target_include_directories(deto PUBLIC ${MPI_CXX_INCLUDE_PATH} lammps/src)
target_link_libraries(deto ${MPI_CXX_LIBRARIES} ${LAMMPS_LIBRARIES} ${JPEG_LIBRARIES} ${PNG_LIBRARIES} ${VTK_LIBRARIES})

if (DETO_WAIT_ATTACH)
  target_compile_definitions(deto PUBLIC DETO_WAIT_ATTACH)
endif()

if (DETO_WAIT_ATTACH_RANK GREATER_EQUAL 0)
  target_compile_definitions(deto PUBLIC DETO_WAIT_ATTACH_RANK=${DETO_WAIT_ATTACH_RANK})
endif()
