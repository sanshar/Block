cmake_minimum_required (VERSION 2.8)
project (Block)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 ")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBLAS -DUSELAPACK -DFAST_MTP -D_HAS_CBLAS -D_HAS_INTEL_MKL -fPIC")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -dynamic")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g ")
if(DEBUG)
  message("Warning: Debug mode")
else(DEBUG)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -funroll-loops -DNDEBUG")
endif(DEBUG)

set(CMAKE_VERBOSE_MAKEFILE OFF)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# set SERIAL on if mpi is not wanted.
if (SERIAL)
  message("Warning: MPI is not used")
else (SERIAL)
  find_package(MPI REQUIRED)
  set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
  include_directories(${MPI_CXX_INCLUDE_PATH})
  link_directories(${MPI_CXX_LIBRARIES})
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_CXX_LINK_FLAGS}")
endif (SERIAL)


#set(BLA_VENDOR "Intel10_64lp_seq")
#Enable Fortran is needed to find the sgemm_ function in some OS.
enable_language(Fortran)
find_package(BLAS)
if (BLA_VENDOR MATCHES "Intel*")
  if (${DISABLE_AVX})
    find_library(BLAS_mkl_def_LIBRARY
      NAMES mkl_def
      PATHS ENV;LD_LIBRARY_PATH)
    set(BLAS_LIBRARIES ${BLAS_LIBRARIES};${BLAS_mkl_def_LIBRARY})
  else()
    find_library(BLAS_mkl_avx_LIBRARY
      NAMES mkl_avx
      PATHS ENV;LD_LIBRARY_PATH)
    set(BLAS_LIBRARIES ${BLAS_LIBRARIES};${BLAS_mkl_avx_LIBRARY})
  endif()
endif()
find_package(LAPACK REQUIRED)

if (SERIAL)
  find_package(Boost REQUIRED COMPONENTS system filesystem serialization)
else(SERIAL)
  find_package(Boost REQUIRED COMPONENTS mpi system filesystem serialization)
  endif(SERIAL)

INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )
link_directories(${Boost_LIBRARY_DIR})
link_libraries( ${Boost_LIBRARIES} )
message("Boost library ${Boost_LIBRARIES}")
        
if(${Boost_VERSION} GREATER "105599")
  add_definitions (-DBOOST_1_56_0)
endif()

add_custom_target(libnewmat.a ALL
                  COMMAND ${CMAKE_MAKE_PROGRAM}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/newmat10
                  COMMENT "Building newmat")
#add_library(${CMAKE_CURRENT_SOURCE_DIR}/newmat10/libnewmat.a STATIC IMPORTED)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/newmat10)
link_directories(${CMAKE_CURRENT_SOURCE_DIR}/newmat10)
link_libraries( ${CMAKE_CURRENT_SOURCE_DIR}/newmat10/libnewmat.a )

include_directories(include)
include_directories(genetic)
include_directories(modules/ResponseTheory)
include_directories(modules/four_index_ops/)
include_directories(modules/generate_blocks/)
include_directories(modules/mps_nevpt/)
include_directories(modules/nevpt2/)
include_directories(modules/npdm/)
include_directories(modules/onepdm/)
include_directories(modules/three_index_ops/)
include_directories(modules/two_index_ops/)
include_directories(modules/twopdm/)
include_directories(molcas)

file(GLOB c_files *.C modules/*/*.C include/*.C molcas/*.C genetic/CrossOver.C genetic/Evaluate.C genetic/GAInput.C genetic/GAOptimize.C genetic/Generation.C genetic/Mutation.C genetic/RandomGenerator.C genetic/ReadIntegral.C )
list(REMOVE_ITEM c_files ${CMAKE_CURRENT_SOURCE_DIR}/main.C)
list(REMOVE_ITEM c_files ${CMAKE_CURRENT_SOURCE_DIR}/OverlapHelement.C)
list(REMOVE_ITEM c_files ${CMAKE_CURRENT_SOURCE_DIR}/Coefficients.C)
list(REMOVE_ITEM c_files ${CMAKE_CURRENT_SOURCE_DIR}/CSFOverlapHelement.C)



add_library(object_files OBJECT ${c_files})
add_executable(block.spin_adapted main.C $<TARGET_OBJECTS:object_files>)
add_executable(OH OverlapHelement.C $<TARGET_OBJECTS:object_files>)
add_executable(COEF Coefficients.C $<TARGET_OBJECTS:object_files>)
add_executable(CSFOH CSFOverlapHelement.C $<TARGET_OBJECTS:object_files>)

target_link_libraries(block.spin_adapted ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries(OH ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries(COEF ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries(CSFOH ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

add_dependencies(object_files libnewmat.a)
#add_dependencies(block.spin_adapted libnewmat.a)
#add_dependencies(OH libnewmat.a)
#add_dependencies(COEF libnewmat.a)
#add_dependencies(CSFOH libnewmat.a)
