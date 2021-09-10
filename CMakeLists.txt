cmake_minimum_required( VERSION 2.5)
if(COMMAND cmake_policy)
    cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/sympiler/cmake)
include(DownloadProject)


set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/lib/)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/lib/)
set(SE_EXTERNAL ${CMAKE_CURRENT_LIST_DIR}/external/)



project(sympiler_eigen)



find_package(METIS OPTIONAL_COMPONENTS)
if(NOT METIS_FOUND)
    set(METIS_INCLUDE_DIR ${LBC_ROOT_DIR}/external/metis/metis_shared-src/include/)
    set(METIS_INCLUDES ${LBC_ROOT_DIR}/external/metis/metis_shared-src/include/)
    set(METIS_LIBRARY ${LBC_ROOT_DIR}/lib/)
    set(METIS_LIBRARIES ${LBC_ROOT_DIR}/lib/)
    set(METIS_LIB metis)
    set(METIS_FOUND TRUE)
else()
    set(METIS_LIB ${METIS_LIBRARIES})
endif()
set(SYMPILER_INC "${CMAKE_CURRENT_SOURCE_DIR}/sympiler/includes/")
set(SYMPILER_LBC_INC "${CMAKE_CURRENT_SOURCE_DIR}/sympiler/lbc/includes/")

set(SYMPILER_LINK_DIR "${CMAKE_CURRENT_SOURCE_DIR}/sympiler/lib/")
set(LBC_LINK_DIR "${CMAKE_CURRENT_SOURCE_DIR}/sympiler/lbc/lib/")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -w -m64")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=x86-64   -mavx2")



download_project(PROJ eigen
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git  # https://github.com/KarypisLab/METIS.git
        PREFIX  ${SE_EXTERNAL}/eigen/
        BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/eigen/
        UPDATE_DISCONNECTED 1
        )

set(EIGEN_INCLUDE_DIR ${eigen_SOURCE_DIR}/)


#find_package(METIS OPTIONAL_COMPONENTS)
link_directories(
    ${SYMPILER_LINK_DIR}
    ${LBC_LINK_DIR}
)

include_directories( ${CMAKE_CURRENT_SOURCE_DIR}includes
        ${SYMPILER_INC}
        ${SYMPILER_LBC_INC}
        ${EIGEN_INCLUDE_DIR}
        )

Set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")

add_subdirectory(sympiler)
add_library(sympiler_eigen
        sympiler_cholesky.cpp)

target_link_libraries(sympiler_eigen
        sym_sparse_blas
        lbc
        ${METIS_LIB}
        )

# Uncomment the following line for building examples
add_subdirectory(example)




