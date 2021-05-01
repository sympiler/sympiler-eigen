cmake_minimum_required(VERSION 2.8)



include_directories(../includes

        ${SYMPILER_SPARSE_BLAS_INC}
        ${SYMPILER_SPARSE_BLASP_INC}
        )

link_directories(
        ${CMAKE_CURRENT_LIST_DIR}/../lib/
)

link_libraries(
        lbc
        sym_sparse_blas
        sym_dense_blas
        sympiler_utils
        ${METIS_LIBRARY} ${METIS_LIBRARIES})

project(sympiler_eigen_demo)
set(SYMPILERD_SOURCE_FILES
        ${METIS_SOURCE_FILES}
        ${ALL_SOURCE_FILES}
        )

add_executable(sym_cholesky_eigen ${SYMPILERD_SOURCE_FILES}
        cholesky_eigen_driver.cpp )
