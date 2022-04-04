#include "bfs.hpp"

#include <vector>
#include <iostream>

// matrix, vertex, girth (current known lowest)
unsigned int bfs(const csr_matrix& m, size_t v, unsigned int g) {
    std::vector<bool> mask(m.size, true); // changed to bool, logic inverted
    std::vector<size_t> vertex_queue(m.size, 0);

    vertex_queue[0] = v;

    auto step_slice_start = vertex_queue.begin();
    auto step_slice_end = step_slice_start + 1;
    auto next_slice = step_slice_end;

    for (unsigned int step = 1; step < g; step++) {
        // bool unchanged = true;
        for (auto i = step_slice_start; i != step_slice_end; i++) {
            const auto row = *i;
            for (auto element_of_row = m.row_ptr[row]; element_of_row < m.row_ptr[row + 1]; element_of_row++) {
                const auto column = m.col_idx[element_of_row];
                if (column == v) return step;//std::cout << "what?" << std::endl;
                if (mask[column]) {
                    // y[column] = x[row];
                    *next_slice = column;
                    mask[column] = false;
                    // unchanged = false;
                    next_slice++;

                    // golfed alternative
                    // mask[*(next_slice++) = column] = false;
                }
            }

            mask[row] = false;
            // x[row] = 0;
        }
        // TODO: check:
        // if (y[v]) return step; // this only handles self-loops, I think. see what this does for 2-cycles
        // if (unchanged) return step;
        // for (int i = 0; i < y.size(); i++) {
        //     if (y[i]) std::cout << i << ",";
        // }
        // std::cout << std::endl;
        step_slice_start = step_slice_end;
        step_slice_end = next_slice;
        // std::swap(x, y);
    }
    return g;
}