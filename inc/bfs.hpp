#ifndef BFS_H
#define BFS_H

#include "csr_matrix.hpp"

// TODO: maybe we can also mask out vertices which we already know the cycle sizes for?

unsigned int bfs(const csr_matrix& m, size_t v, unsigned int g);

#endif