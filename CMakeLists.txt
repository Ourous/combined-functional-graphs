cmake_minimum_required(VERSION 3.0.0)
project(girth-of-quadratic-graphs VERSION 1.0.0)

set(CMAKE_CXX_FLAGS "-std=c++2a -Wno-newline-eof -Wno-ctad-maybe-unsupported -Wno-c++98-compat -Wno-c++98-compat-pedantic -Werror -Wshadow -Wconversion -Wall -Wextra -pthread")

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3 -march=native")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -D_EG_DEBUG -Og -pg")

include_directories(${CMAKE_SOURCE_DIR}/inc)

find_library(GMP_LIBRARY gmp REQUIRED)

add_subdirectory(lib)
add_subdirectory(src)

add_executable(girthtool main.cpp)
target_link_libraries(girthtool bfs thread_pool gmp)