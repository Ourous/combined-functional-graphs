#ifndef BFS_H
#define BFS_H

#include "csr_matrix.hpp"

// TODO: can also mask out vertices which we already know the cycle sizes for?

int bfs(const csr_matrix& m, size_t v, int g);

#endif