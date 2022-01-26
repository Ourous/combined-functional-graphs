#include <iostream>
#include <chrono>

#include <csr_matrix.hpp>
#include <bfs.hpp>


int main(int, char**) {
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

    std::set<std::pair<int, int>> points;

    uintmax_t p = 7919;

    for (auto x = 0; x < p; x++) {
        points.insert(std::make_pair((x * x + 6) % p, x));
    }
    csr_matrix m(p, points);
    int g = INT32_MAX;
    for (int v = 0; v < p; v++) {
        g = std::min(g, bfs(m, v, g));
        // auto girth = bfs(m, v, g);
        // std::cout << "from " << v << ": " << girth << std::endl;
    }

    std::cout << g << std::endl;

    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    std::cerr << "Finished in: " << duration << " ms" << std::endl;
}
