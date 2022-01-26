#include "csr_matrix.hpp"

csr_matrix::csr_matrix(size_t size, std::set<std::pair<size_t, size_t>> coords) : size(size), col_idx(), row_ptr() {
    col_idx.reserve(coords.size());
    row_ptr.reserve(size + 1);
    for (size_t y = 0; y < size; y++) {
        bool first_nz = true;
        for (size_t x = 0; x < size; x++) { // TODO: c++20 -> set.contains
            if (coords.count(std::make_pair(x, y))) {
                if (first_nz) { // TODO: unwrap this into two partial loops instead
                    row_ptr.push_back(col_idx.size());
                    first_nz = false;
                }
                col_idx.push_back(x);
            }
        }
        if (first_nz) {
            size_t back = row_ptr.back();
            row_ptr.push_back(back);
        }
    }
    row_ptr.push_back(col_idx.size());
}