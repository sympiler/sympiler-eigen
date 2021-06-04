cmake_minimum_required(VERSION 2.8)



include_directories(../includes
        )

link_directories(
        ${CMAKE_CURRENT_LIST_DIR}/../lib/
)


project(sympiler_eigen_demo)
set(SYMPILERD_SOURCE_FILES
        ${METIS_SOURCE_FILES}
        ${ALL_SOURCE_FILES}
        )
add_executable(sympiler_eigen_demo ${SYMPILERD_SOURCE_FILES}
        cholesky_eigen_driver.cpp )
target_link_libraries(sympiler_eigen_demo
        sympiler_eigen
        metis
)
