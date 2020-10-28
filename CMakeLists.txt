# CMakeLists.txt for event package. It creates a library with dictionary and a main program
cmake_minimum_required(VERSION 2.8.10)
project (testMapPointing)
set(libname "TestMP")

set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_COLOR_MAKEFILE ON)

# If ROOTSYS is defined in the usual way then this should be sufficient to find the ROOT CMake things
if(DEFINED $ENV{ROOTSYS})
  list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
# If not then we can access the variables the ROOT CMake stuff directly, fingers crossed
else()
  execute_process(COMMAND root-config --etcdir OUTPUT_VARIABLE ROOT_ETC_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${ROOT_ETC_DIR}/cmake)
endif()

if(DEFINED ENV{ANITA_UTIL_INSTALL_DIR})
    message("ANITA_UTIL_INSTALL_DIR is set to $ENV{ANITA_UTIL_INSTALL_DIR}")
    set(UTIL_LIB_DIR $ENV{ANITA_UTIL_INSTALL_DIR}/lib)
    set(UTIL_INC_DIR $ENV{ANITA_UTIL_INSTALL_DIR}/include)
    set(UTIL_BIN_DIR $ENV{ANITA_UTIL_INSTALL_DIR}/bin)
    set(UTIL_SHARE_DIR $ENV{ANITA_UTIL_INSTALL_DIR}/share)
    set(LD_UTIL $ENV{ANITA_UTIL_INSTALL_DIR}/lib)
    set(INC_UTIL $ENV{ANITA_UTIL_INSTALL_DIR}/include)
elseif(DEFINED ENV{ARA_UTIL_INSTALL_DIR})
    message("ARA_UTIL_INSTALL_DIR is set to $ENV{ANITA_UTIL_INSTALL_DIR}")
    set(UTIL_LIB_DIR $ENV{ARA_UTIL_INSTALL_DIR}/lib)
    set(UTIL_INC_DIR $ENV{ARA_UTIL_INSTALL_DIR}/include)
    set(UTIL_BIN_DIR $ENV{ARA_UTIL_INSTALL_DIR}/bin)
    set(UTIL_SHARE_DIR $ENV{ARA_UTIL_INSTALL_DIR}/share)
    set(LD_UTIL $ENV{ARA_UTIL_INSTALL_DIR}/lib)
    set(INC_UTIL $ENV{ARA_UTIL_INSTALL_DIR}/include)
else()
    message("Will install to /usr/local ANITA_UTIL_INSTALL_DIR not set")
    set(UTIL_LIB_DIR /usr/local/lib)
    set(UTIL_INC_DIR /usr/local/include)
    set(UTIL_BIN_DIR /usr/local/bin)
    set(UTIL_SHARE_DIR /usr/local/share)
    set(LD_UTIL /usr/local/lib)
    set(INC_UTIL /usr/local/include)
endif()

option( FFTTOOLS_ENABLE_OPENMP "If FFTtools was compiled with OpenMP support, enable this option to allow multiple FFT's to be computed in parallel") 
if (FFTTOOLS_ENABLE_OPENMP)
  add_definitions(-DFFTTOOLS_COMPILED_WITH_OPENMP)
endif()

option ( ENABLE_OPENMP " Enable OpenMP for multicore usage") 
if(ENABLE_OPENMP) 
  FIND_PACKAGE( OpenMP REQUIRED)
  if(OPENMP_FOUND)
    message("OPENMP FOUND")
    add_definitions (-DUCORRELATOR_OPENMP) 
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  else() 
    message(FATAL_ERROR "you tried to use openmp but Cmake couldn't find compiler support") 
  endif() 
endif()

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS MathMore MathCore RIO Hist Tree Net Minuit Minuit2 TreePlayer)
message("ROOT_INCLUDE_DIRS is set to ${ROOT_INCLUDE_DIRS}")
#find_package(FFTW REQUIRED)
find_package(ZLIB REQUIRED)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)

message("ROOT_VERSION is set to ${ROOT_VERSION}")
#if( ${ROOT_VERSION} VERSION_GREATER "5.99/99")
#    message("Using ROOT_VERSION 6")
#    include(${ROOT_USE_FILE})
#else()
##    include(RootNewMacros)
#     add_definitions(${ROOT_DEFINITIONS})
#endif()

message("ROOT_INCLUDE_DIRS is set to ${ROOT_INCLUDE_DIRS}")

include_directories(${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR}/include ${FFTW_INCLUDES} ${ZLIB_INCLUDES} ${UTIL_INC_DIR})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  ${ROOT_CXX_FLAGS}")
add_definitions(${ROOT_CXX_FLAGS} -DSummarySelector_cxx)
link_directories(${UTIL_LIB_DIR})

#file(GLOB HEADER_FILES ${PROJECT_SOURCE_DIR}/include/*.h)
#file(GLOB_RECURSE SOURCE_FILES ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB HEADER_FILES "include/*.h")
file(GLOB_RECURSE SOURCE_FILES src/*.cc)

#set(DICTNAME G__${libname})

#ROOT_GENERATE_DICTIONARY(${DICTNAME} ${HEADER_FILES} LINKDEF LinkDef.h)

#---Create a shared library with generated dictionary
add_library(${libname} SHARED ${SOURCE_FILES})
target_link_libraries(${libname} ${ZLIB_LIBRARIES} AnitaAnalysis AnitaEvent RootFftwWrapper UCorrelator ${ROOT_LIBRARIES} ${FFTW_LIBRARIES})
#add_library(${libname} SHARED ${SOURCE_FILES} ${DICTNAME}.cxx)
#target_link_libraries(${libname} ${ZLIB_LIBRARIES} AnitaAnalysis AnitaEvent AnitaCorrelator RootFftwWrapper UCorrelator ${ROOT_LIBRARIES} ${FFTW_LIBRARIES})

#if( ${ROOT_VERSION} VERSION_GREATER "5.99.99")
#    add_custom_target(${DICTNAME}.pcm DEPENDS ${DICTNAME})
#endif()

message("UTIL_INC_DIR is set to ${UTIL_INC_DIR}")
#message("DICTNAME is set to ${DICTNAME}" )
message("PROJECT_BINARY_DIR is set to ${PROJECT_BINARY_DIR}")


set(ANITA_LIBS AnitaEvent AnitaAnalysis AnitaAnalysisTools AnitaCorrelator RootFftwWrapper UCorrelator)

set(BINARIES makeDeltaTMapFiles makeAntSphericalCosineProductsFiles testMapPointing)

#add_library(testMapPointingLib SHARED ${SOURCE_FILES})
#target_link_libraries(testMapPointingLib ${SOURCE_FILES} ${ZLIB_LIBRARIES} ${ROOT_LIBRARIES} ${FFTW_LIBRARIES} ${ANITA_LIBS})

FOREACH(binary ${BINARIES})
  MESSAGE(STATUS "Process file: ${binary}")
  add_executable(${binary} ${binary}.cxx)
  target_link_libraries(${binary} ${libname})
#  target_link_libraries(${binary} ${ZLIB_LIBRARIES} ${ANITA_LIBS} ${ROOT_LIBRARIES} ${FFTW_LIBRARIES})
ENDFOREACH(binary)