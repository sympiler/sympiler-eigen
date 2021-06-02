cmake_minimum_required( VERSION 2.5)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/sympiler/cmake)
include(DownloadProject)


set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/lib/)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/lib/)
project(sympiler_eigen)


set(SYMPILER_UTIL_INC "${CMAKE_CURRENT_SOURCE_DIR}/sympiler/utils/includes")
set(SYMPILER_UTILE_INC "${CMAKE_CURRENT_SOURCE_DIR}/sympiler/utils/external/includes")
set(SYMPILER_SPARSE_BLAS_INC "${CMAKE_CURRENT_SOURCE_DIR}/sympiler/sparse_blas/includes")
set(SYMPILER_SPARSE_BLASP_INC "${CMAKE_CURRENT_SOURCE_DIR}/sympiler/sparse_blas/includes/parsy")
set(SYMPILER_LBC_INC "${CMAKE_CURRENT_SOURCE_DIR}/sympiler/lbc/src/includes")

set(SYMPILER_LINK_DIR "${CMAKE_CURRENT_SOURCE_DIR}/sympiler/lib/")
set(LBC_LINK_DIR "${CMAKE_CURRENT_SOURCE_DIR}/sympiler/lbc/lib/")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -w -m64")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=x86-64   -mavx2")


#find_package(METIS OPTIONAL_COMPONENTS)
link_directories(
    ${SYMPILER_LINK_DIR}
        ${LBC_LINK_DIR}
)

include_directories( ${CMAKE_CURRENT_SOURCE_DIR}includes
        ${SYMPILER_SPARSE_BLAS_INC}
        ${SYMPILER_SPARSE_BLASP_INC}
        )
Set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")

add_subdirectory(sympiler)
add_library(sympiler_eigen
        sympiler_cholesky.cpp)
target_link_libraries(sympiler_eigen
        sym_sparse_blas
        sym_dense_blas
        lbc
        sympiler_utils
        )
# Uncomment the following line for building examples
add_subdirectory(example)




