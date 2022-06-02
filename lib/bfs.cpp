#include "bfs.hpp"

#include <vector>

unsigned int bfs(const csr_matrix& m, size_t v, unsigned int g) {
    std::vector<bool> mask(m.size, true);
    std::vector<size_t> vertex_queue(m.size, 0);

    vertex_queue[0] = v;

    auto step_slice_start = vertex_queue.begin();
    auto step_slice_end = step_slice_start + 1;
    auto next_slice = step_slice_end;

    for (unsigned int step = 1; step < g; step++) {

        for (auto i = step_slice_start; i != step_slice_end; i++) {

            const auto row = *i;

            for (auto element_of_row = m.row_ptr[row]; element_of_row < m.row_ptr[row + 1]; element_of_row++) {

                const auto column = m.col_idx[element_of_row];

                if (column == v) return step;

                if (mask[column]) {
                    *next_slice = column;
                    mask[column] = false;
                    next_slice++;
                }
            }

            mask[row] = false;
        }

        step_slice_start = step_slice_end;
        step_slice_end = next_slice;
    }

    return g;
}