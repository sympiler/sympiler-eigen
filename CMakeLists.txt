cmake_minimum_required( VERSION 2.5)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/sympiler/cmake)
include(DownloadProject)


project(sympiler_eigen)


set(SYMPILER_UTIL_INC "${CMAKE_CURRENT_SOURCE_DIR}/sympiler/utils/includes")
set(SYMPILER_UTILE_INC "${CMAKE_CURRENT_SOURCE_DIR}/sympiler/utils/external/includes")
set(SYMPILER_SPARSE_BLAS_INC "${CMAKE_CURRENT_SOURCE_DIR}/sympiler/sparse_blas/includes")
set(SYMPILER_SPARSE_BLASP_INC "${CMAKE_CURRENT_SOURCE_DIR}/sympiler/sparse_blas/includes/parsy")
set(SYMPILER_LBC_INC "${CMAKE_CURRENT_SOURCE_DIR}/sympiler/lbc/src/includes")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -w -m64")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=x86-64   -mavx2")

add_subdirectory(sympiler)
add_subdirectory(example)




