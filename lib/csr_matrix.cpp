#include "csr_matrix.hpp"

csr_matrix::csr_matrix(size_t size, const std::map<size_t, std::set<size_t>>& coords) : size(size), col_idx(), row_ptr() {
    row_ptr.reserve(size + 1);

    for (size_t y = 0; y < size; y++) {

        if (coords.contains(y)) {

            row_ptr.push_back(col_idx.size());
            col_idx.insert(col_idx.end(), coords.at(y).begin(), coords.at(y).end());

        }
        else row_ptr.push_back(row_ptr.size() ? row_ptr.back() : 0);
    }

    row_ptr.push_back(col_idx.size());
}