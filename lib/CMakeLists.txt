include_directories(${CMAKE_SOURCE_DIR}/inc)

add_library(thread_pool thread_pool.cpp)

add_library(quadratic quadratic.cpp)

add_library(csr_matrix csr_matrix.cpp)

add_library(bfs csr_matrix bfs.cpp)
