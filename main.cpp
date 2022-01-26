#include <iostream>
#include <chrono>

#include <vector>
#include <set>
#include <utility>
#include <map>
#include <bitset>
#include <algorithm>

#include <csr_matrix.hpp>
#include <bfs.hpp>
#include <pari.h>
#include <gmp.h>

template<typename T>
std::set<T> set_merge_linear_time(const std::set<T>& s1, const std::set<T>& s2) {
    std::vector<T> v;
    std::merge(s1.begin(), s1.end(), s2.begin(), s2.end(), back_inserter(v));
    // TODO: it's possible to do this while checking every time for girth = 2, in still linear time
    // .. because set access is constant-amortized
    return std::set<T>(v.begin(), v.end());
}


// int determine_girth_of(const std::set<size_t>) {
//     if (size_t)
// }

void compute_girth_stats(size_t p) {
    // initialize the various function sets
    std::map<size_t, std::set<std::pair<size_t, size_t>>> quadratics;
    std::vector<size_t> others;



    std::set<size_t> ones;
    // TODO: store single twos in different set to speed up subset checking
    std::set<std::set<size_t>> twos; // never need to store anything if a subset is already in here
    // std::map<std::set<size_t>, size_t> others; // once we fill out 3s we can delete 2s and so-on
    // std::map<std::set<size_t>, int> girths;
    // std::map<int,std::set<std::set<size_t>>> girths;
    std::vector<std::vector<std::vector<size_t>>> girths(p, std::vector<std::vector<size_t>>());

    for (size_t a = 0; a < p; a++) {
        std::set<std::pair<size_t, size_t>> points;
        bool girth_is_one = false, girth_maybe_two = false;
        for (size_t x = 0; x < p && !girth_is_one; x++) { // works for p < 2^63
            size_t y = (x * x + a) % p;
            if (x == y) girth_is_one = true;
            else if (!girth_maybe_two) {
                if (x == (y * y + a) % p) girth_maybe_two = true;
                else
                    points.insert(std::make_pair(y, x));
            }
        }
        if (girth_is_one) {
            ones.insert(a);
            // girths.insert({ { a }, 1 });
            // girths[1].push_back({ a });
        }
        else if (girth_maybe_two) {
            twos.insert({ a });
            // girths.insert({ { a }, 2 });
            // girths[2].push_back({ a });
        }
        else {
            quadratics.insert({ a, points });
            others.push_back(a);
        }

    }


    for (auto [a, points] : quadratics) {
        // for (size_t a = 0; a < p; a++) {
            // if (ones.count(a)) {
            //     std::cout << a << ": " << 1 << std::endl;
            //     continue;
            // }
            // if (twos.count({ a })) {
            //     std::cout << a << ": " << 2 << std::endl;
            //     continue;
            // }
        csr_matrix m(p, points);
        int g = p;
        for (size_t v = 0; v < p; v++) {
            g = std::min(g, bfs(m, v, g));
            // auto girth = bfs(m, v, g);
            // std::cout << "from " << v << ": " << girth << std::endl;
        }
        // girths.insert({ { a }, g });
        girths[g].push_back({ a });

        std::cout << a << ": " << g << std::endl;
    }


    for (int s = 2; s <= quadratics.size(); s++) {
        std::vector<bool> combination(quadratics.size(), false);
        std::fill(combination.begin(), combination.begin() + s, true);
        do {
            std::vector<size_t> coefficients;
            for (int i = 0; i < combination.size(); i++) {
                if (!combination[i]) continue;
                coefficients.push_back(others[i]);
            }
            int g = 0;
            [&] {
                for (; g < girths.size(); g++) {
                    for (auto subset : girths[g]) {
                        if (std::includes(coefficients.begin(), coefficients.end(), subset.begin(), subset.end())) {
                            return;
                        }
                    }
                }
            }();
            if (g <= s) continue;
            if (g > s) {
                std::cout << "g: " << g << ", s: " << s << std::endl;
                std::set<std::pair<size_t, size_t>> combined_function;
                for (auto a : coefficients) {
                    // combined_function = set_merge_linear_time(quadratics[a], combined_function);
                    combined_function.insert(quadratics[a].begin(), quadratics[a].end());
                }

                // for g = 1, there is a subset with (s=1,g=1)
                // for g = 2, there is a subset with (s<2,g=2)
                // for g = 3, there is a subset with (s<3,g=3)
                // for g = g`, there is a subset with (s<g`,g=g`)

                // BIG IDEA: 

                // for s = 3, if g < 3 we will know via subset
                // for s = 4, if g < 4 we will know via subset
                // for s = 5, if g < 5 we will know via subset

                // once we retreive the initial girth via checking subsets
                // if the initial g = s: that is the final girth

                // so we only store subsets where g >= s

                // I *think* we don't need to check where g == s and we can assume it remains the same
                // TODO: CHECK THE ABOVE LINE; IMPLEMENT IF TRUE

                // SECOND IDEA:

                // we can remove single functions with g = N when s > N
                // TODO: TEST WITH DATA ^

                csr_matrix m(p, combined_function);
                int initial_g = g;
                for (size_t v = 0; v < p; v++) g = std::min(g, bfs(m, v, g));
                if (initial_g > g) girths[g].push_back(coefficients);
            }
            // girths.insert({ coefficients, g });
            for (auto i = coefficients.begin(); i != coefficients.end();) {
                std::cout << *i;
                if (i++ == coefficients.end()) break;
                else std::cout << ",";
            }
            std::cout << ": " << g << std::endl;

        } while (std::next_permutation(combination.rbegin(), combination.rend()));
    }
}

int main(int, char**) {
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

    size_t p = 101;
    //97 = 17 seconds with twos, 5 seconds with full lookup
    //

    compute_girth_stats(p);

    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    std::cerr << "Finished in: " << duration << " ms" << std::endl;
}
