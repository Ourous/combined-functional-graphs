#ifndef CSR_MATRIX_H
#define CSR_MATRIX_H
// think on this
// we don't need to store the entries, they're all =1
// the col_index is just a list of essentially the (output) of the function
// the row_pointer is 0s, 1s, 2s, 3s ... (p)s

#include <vector>
#include <set>
#include <utility>

struct csr_matrix {
    int size;
    std::vector<int> col_idx;
    std::vector<int> row_ptr;
    csr_matrix(int, std::set<std::pair<int, int>>);
};

#endif