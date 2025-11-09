cmake_minimum_required (VERSION 3.0...3.5)                                            
project(mugast CXX)                                                              
set(CMAKE_BUILD_TYPE Release)                                                   
include(CMakePackageConfigHelpers)

set(CMAKE_BINARY_OUTPUT_DIRECTORY  ${CMAKE_INSTALL_PREFIX}/bin )
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_INSTALL_PREFIX}/lib )
set(CMAKE_INCLUDE_OUTPUT_DIRECTORY ${CMAKE_INSTALL_PREFIX}/include )
set(CMAKE_INSTALL_LIBDIR ${CMAKE_LIBRARY_OUTPUT_DIRECTORY} )
### nptool ###
find_package(nptool QUIET)

include_directories($ENV{NPTOOL}/include)
include_directories($ENV{NPTOOL_HOME}/$ENV{NPTOOL_ENV}/include)

link_directories($ENV{NPTOOL}/lib)
link_directories(../install/lib)
link_directories(../../install/lib)
link_directories($ENV{NPTOOL_HOME}/$ENV{NPTOOL_ENV}/lib)

if(DEFINED ENV{MFMLIB})
  SET(MFMLIB $ENV{MFMLIB})
  SET(MFM_DIR $ENV{MFMLIB})
  MESSAGE(STATUS "The libraries and software will be build using MFM libraries")
  SET(MFM_LIBRARY_DIR ${MFM_DIR}/lib)
  # Set MFM_BIN_DIR
  SET(MFM_BIN_DIR ${MFM_DIR}/bin)
  # Set MFM_INCLUDE_DIR
  SET(MFM_INCLUDES ${MFM_DIR}/include)
  SET(MFM_INCLUDE_DIRS  ${MFM_DIR}/include )
  SET(MFM_LIBRARIES "-L${MFM_LIBRARY_DIR} -lMFM")
  LIST(APPEND COMPDEF "WITH_MFM")
  LIST(APPEND MFMCOMPLIB ${MFM_LIBRARIES})
  LIST(APPEND INCLUDE_DIR ${MFM_INCLUDE_DIRS})
  if(EXISTS "${MFM_DIR}/lib/libMFM${CMAKE_SHARED_LIBRARY_SUFFIX}")
    set(MFM_FOUND true)
  else()
    set(MFM_FOUND false)
  endif()
endif()

### ROOT ###
# find package 
find_package(ROOT QUIET)
if(${ROOT_FOUND})
  message(STATUS "ROOT found: ${ROOT_VERSION}" )
  # use this definition to isolate ROOT dependent code
  add_compile_definitions(ROOT_FOUND)
  # add usefull root fonction for dictionnary definition
  include(${ROOT_USE_FILE})
  # add include directory
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
else()
  message(STATUS "ROOT not found: compiling with no ROOT support" )
endif()

### GEANT4 ###
find_package(Geant4 QUIET)
if(${Geant4_FOUND})
  message(STATUS "Geant4 found: ${Geant4_VERSION}" )
  # use this definition to isolate Geant4 dependent code
  add_compile_definitions(Geant4_FOUND)
  # add usefull root fonction for dictionnary definition
  include(${Geant4_USE_FILE})
  # add include directory
  include_directories(${Geant4_INCLUDE_DIRS})
else()
  message(STATUS "Geant4 not found: compiling with no Geant4 support" )
endif()

set (duplicate_include_dir ${CMAKE_CURRENT_SOURCE_DIR}/src PARENT_SCOPE)
include_directories("src")
add_subdirectory("src")
