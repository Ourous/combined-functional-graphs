#ifndef CSR_MATRIX_H
#define CSR_MATRIX_H
// think on this
// we don't need to store the entries, they're all =1
// the col_index is just a list of essentially the (output) of the function
// the row_pointer is 0s, 1s, 2s, 3s ... (p)s

#include <vector>
#include <set>
#include <utility>
#include <map>

struct csr_matrix {
    size_t size;
    std::vector<size_t> col_idx;
    std::vector<size_t> row_ptr;
    csr_matrix(size_t, std::set<std::pair<size_t, size_t>>); // x,y
    csr_matrix(size_t, const std::map<size_t, std::set<size_t>>&); // y,x
};

#endif