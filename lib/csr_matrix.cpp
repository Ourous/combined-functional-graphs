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
            row_ptr.push_back(row_ptr.size() ? row_ptr.back() : 0);
        }
    }
    row_ptr.push_back(col_idx.size());
}

// TODO: see if I can do with with: const std::unordered_map<size_t, std::vector<size_t>>& coords
csr_matrix::csr_matrix(size_t size, std::map<size_t, std::set<size_t>> coords) : size(size), col_idx(), row_ptr() {
    //assuimes all sets are non-empty
    row_ptr.reserve(size + 1);
    for (size_t y = 0; y < size; y++) {
        if (coords.contains(y)) {
            row_ptr.push_back(col_idx.size());
            // for (size_t x : coords[y]) col_idx.push_back(x); // TODO optimize this
            col_idx.insert(col_idx.end(), coords[y].begin(), coords[y].end());
        }
        else row_ptr.push_back(row_ptr.size() ? row_ptr.back() : 0);
    }
    row_ptr.push_back(col_idx.size());
}