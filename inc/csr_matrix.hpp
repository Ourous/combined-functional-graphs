#ifndef CSR_MATRIX_H
#define CSR_MATRIX_H

#include <vector>
#include <set>
#include <utility>
#include <map>

struct csr_matrix {
    size_t size;
    std::vector<size_t> col_idx;
    std::vector<size_t> row_ptr;
    csr_matrix(size_t, const std::map<size_t, std::set<size_t>>&); // y,x
};

#endif