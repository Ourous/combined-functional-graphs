include_directories(${CMAKE_SOURCE_DIR}/inc)

add_library(quadratic quadratic.cpp)

add_library(csr_matrix csr_matrix.cpp)

add_library(bfs csr_matrix bfs.cpp)
