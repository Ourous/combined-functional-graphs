cmake_minimum_required(VERSION 3.0.0)
project(girth-of-quadratic-graphs VERSION 0.1.0)

set(CMAKE_CXX_FLAGS "-std=c++2a -Wno-newline-eof -Wno-ctad-maybe-unsupported -Wno-unused-command-line-argument -Wno-c++98-compat -Wno-c++98-compat-pedantic -Werror -Wshadow -Wconversion -Wall -Wextra -pthread -lgmp -lgmpxx")

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -Ofast")
# set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3 -fgcse-sm -fgcse-las")
# set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -Wc++11-compat -D_EG_DEBUG -pg")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -D_EG_DEBUG -Og -pg")

include(CTest)
enable_testing()

include_directories(${CMAKE_SOURCE_DIR}/inc)

include_directories(/usr/include/x86_64-linux-gnu/pari)

find_library(PARI_LIBRARY pari REQUIRED)
find_library(GMP_LIBRARY gmp REQUIRED)

add_subdirectory(lib)
add_subdirectory(src)

add_executable(girthtool main.cpp)
target_link_libraries(girthtool bfs thread_pool gmp)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
