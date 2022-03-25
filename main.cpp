#include <iostream>
#include <chrono>

#include <thread>

#include <vector>
#include <set>
#include <utility>
#include <map>
#include <bitset>
#include <algorithm>
#include <shared_mutex>
#include <functional>
#include <queue>
#include <cassert>

#include <csr_matrix.hpp>
#include <bfs.hpp>
#include <gmpxx.h>
#include <thread_pool.hpp>

template<typename T>
std::set<T> set_merge_linear_time(const std::set<T>& s1, const std::set<T>& s2) {
    std::vector<T> v;
    std::merge(s1.begin(), s1.end(), s2.begin(), s2.end(), back_inserter(v));
    // TODO: it's possible to do this while checking every time for girth = 2, in still linear time
    // .. because set access is constant-amortized
    return std::set<T>(v.begin(), v.end());
}


uintmax_t binomial(uintmax_t n, uintmax_t k) {
    if (n == k) return 1;
    if (n < k) return 0;
    mpz_class out;
    mpz_bin_uiui(out.get_mpz_t(), n, k);
    if (!out.fits_ulong_p()) std::cout << "no" << std::endl;
    // segfault either because of GMP or because of combination iteration overrun (off-by-one error)
    return out.get_ui();
}

std::vector<size_t> find_k_combination(size_t k, size_t position) {
    std::vector<size_t> combination_indices;
    combination_indices.reserve(k);
    for (; k > 0; k--) {
        size_t n = k - 1;
        while (binomial(n + 1, k) <= position) n++;
        position -= binomial(n, k);
        combination_indices.push_back(n);
    }
    return combination_indices;
}

void compute_girth_stats(size_t p) {
    // initialize the various function sets

    std::map<size_t, std::map<size_t, std::set<size_t>>> quadratics;

    // quadratic_coefficients[index] = coefficient (in order)
    std::vector<size_t> quadratic_coefficients;

    std::vector<size_t> singles(p, 0); // a -> g
    // girth-ordered sets of subsets that dictate a girth
    std::vector<std::vector<std::vector<size_t>>> girths(p, std::vector<std::vector<size_t>>());

    for (size_t a = 0; a < p; a++) {
        std::map<size_t, std::set<size_t>> points;
        bool girth_is_one = false, girth_maybe_two = false;
        for (size_t x = 0; x < p && !girth_is_one; x++) { // works for p < 2^63
            size_t y = (x * x + a) % p;
            if (x == y) girth_is_one = true;
            else if (!girth_maybe_two) {
                if (x == (y * y + a) % p) girth_maybe_two = true;
                else {
                    if (points.contains(x)) points[x].insert(y);
                    else points[x] = { y };
                }
            }
        }
        if (girth_is_one) {
            singles[a] = 1;
        }
        else if (girth_maybe_two) {
            singles[a] = 2;
        }
        else {
            quadratics.insert({ a, points });
            quadratic_coefficients.push_back(a);
        }
    }

    unsigned int limit_g = 0;

    for (auto [a, points] : quadratics) {
        csr_matrix m(p, points);
        unsigned int g = p;
        for (size_t v = 0; v < p; v++) {
            g = std::min(g, bfs(m, v, g));
        }
        limit_g = std::max(g, limit_g);
        singles[a] = g;
        girths[g].push_back({ a });

        // std::cout << a << ": " << g << std::endl;
    }

    girths.resize(limit_g + 1);

    std::mutex girths_mutex;
    size_t t = std::thread::hardware_concurrency();
    thread_pool pool(t);
    std::vector<std::queue<std::pair<unsigned int, std::vector<size_t>>>> coefficient_queues_2nd_pass(t, std::queue<std::pair<unsigned int, std::vector<size_t>>>());
    std::vector<std::mutex> thread_output_locks(t);
    std::vector<std::condition_variable> thread_notifiers(t);
    for (size_t s = 2; s <= quadratics.size(); s++) {

        // todo: split cofficients into (t) queues that can be taken out of by the (t) threads
        //.. and either discarded by the threads or stored in another queue
        //.. when all coefficients are initially examined against the girth subset lookup
        //.. (possibly redistribute the members of these queues evenly first)
        //.. stop these (t) threads and start another (t) threads on the (t) new queues
        //.. which will perform the girth computation and store it until they're done all the work in their queue
        //.. after which they all write their results into the shared girths

        // start 1st pass threads here
        unsigned int max_g = 0;
        std::vector<size_t> max_g_v;

        auto total_combinations = binomial(quadratics.size(), s);
        // std::cout << " TOTAL : " << total_combinations << std::endl;
        auto combinations_per_job = total_combinations / t;
        auto leftover_combinations = total_combinations % t;
        // std::cout << " PER: " << combinations_per_job << std::endl;
        size_t start_combination_id = 0;

        for (size_t which_t = 0; which_t < t; which_t++) {

            auto combinations_in_job = combinations_per_job + ((leftover_combinations > which_t) ? 1 : 0);

            if (combinations_in_job > 0) pool.do_job(
                // I think the order of variable declaration and memory layout in this function is *really* important
                //.. and should be analyzed further: moving the declaration of coefficients around seems to change the runtime by a whole second for 601
                [&, which_t, start_combination_id, combinations_in_job] {

                    unsigned int local_max_g = 0;
                    std::vector<size_t> local_max_g_v;

                    std::vector<bool> combination_indices(quadratics.size(), false);
                    // std::fill(combination.begin(), combination.begin() + s, true);
                    for (auto i : find_k_combination(s, start_combination_id)) { // TODO: just use k-combination index enumeration instead, it's faster
                        combination_indices[i] = true;
                    }

                    std::vector<size_t> coefficients(s);

                    for (unsigned long i = 0; i < combinations_in_job; i++, std::next_permutation(combination_indices.rbegin(), combination_indices.rend())) {
                        size_t coefficient_index = 0;
                        for (size_t i = 0; i < combination_indices.size(); i++) {
                            if (!combination_indices[i]) continue;
                            // coefficients.push_back(quadratic_coefficients[i]);
                            coefficients[coefficient_index++] = quadratic_coefficients[i];
                        }

                        unsigned int g = 0; // memory locality of these subsets may have a massive performance impact
                        [&] {
                            for (; g < girths.size(); g++) {
                                // if we order the subsets lexicographically we can perform a binary search for the range in which subsets that could be inside
                                //.. our coefficients can exist in with reference to the first and last elements of our coefficients
                                //.. eg we can consider only the range [subset.begin() >= coefficients.begin(), subset.end() <= coefficients.end()]
                                for (const auto& subset : girths[g]) {
                                    // TODO: check if std::includes starts with a binary search for the elements at the ends of the ranges
                                    // maybe if I sort `subset` lexicographically?
                                    if (std::includes(coefficients.begin(), coefficients.end(), subset.begin(), subset.end())) { // there is an idea that if I find every possible combination with a coefficient already, I can discard the coefficient
                                        return;
                                    }
                                }
                            }
                        }();

                        if (g <= s) {
                            if (g > local_max_g) {
                                local_max_g = g;
                                local_max_g_v = coefficients;
                            }
                        }
                        else {
                            auto lock = std::scoped_lock(thread_output_locks[which_t]);
                            coefficient_queues_2nd_pass[which_t].emplace(std::pair{ g, coefficients });
                            // maybe also store g here
                        }
                    }

                    {
                        auto lock = std::scoped_lock(girths_mutex);
                        if (local_max_g > max_g) {
                            max_g = local_max_g;
                            max_g_v = local_max_g_v;
                        }
                    }
                }
            );
            start_combination_id += combinations_in_job;
        }


        pool.wait();
        // for (auto& q : coefficient_queues_2nd_pass) {
        //     std::cout << q.size() << std::endl;
        // }

        // do second pass here
        // redistribute queues

        for (size_t which_t = 0; which_t < t; which_t++) {
            if (coefficient_queues_2nd_pass[which_t].empty()) continue;
            pool.do_job(
                [&, which_t] {
                    std::vector<size_t> coefficients;
                    unsigned int g;
                    while (true) {
                        {
                            // maybe use cond_variable here, or scoped_lock
                            if (coefficient_queues_2nd_pass[which_t].empty()) return;

                            auto arguments = std::move(coefficient_queues_2nd_pass[which_t].front());
                            coefficient_queues_2nd_pass[which_t].pop();
                            coefficients = std::move(arguments.second);
                            g = arguments.first;
                        }


                        std::map<size_t, std::set<size_t>> combined_function;
                        for (auto a : coefficients) {
                            // TODO: use STL stuff for speed?
                            for (const auto& [x, y] : quadratics[a]) {
                                if (combined_function.contains(x)) combined_function[x].insert(y.begin(), y.end());
                                else combined_function.emplace(std::pair{ x, y });
                            }
                        }

                        csr_matrix m(p, std::move(combined_function));
                        unsigned int initial_g = g;
                        // initially, we can just make this paralell
                        //.. give each call to BFS a write address that isn't contended or something
                        for (size_t v = 0; v < p && g > s; v++) g = std::min(g, bfs(m, v, g));
                        if (initial_g > g) { // locking on both branches does something funky to the performance
                            {
                                auto lock = std::scoped_lock(girths_mutex);
                                girths[g].push_back(coefficients);
                                if (g > max_g) {
                                    max_g = g;
                                    max_g_v = coefficients;
                                }
                            }
                            // for (auto i = coefficients.begin(); i != coefficients.end();) {
                            //     std::cout << *i;
                            //     if (i++ == coefficients.end()) break;
                            //     else std::cout << ",";
                            // }
                            // std::cout << ": " << g << std::endl;
                        }
                        else {
                            auto lock = std::scoped_lock(girths_mutex);
                            if (g > max_g) {
                                max_g = g;
                                max_g_v = coefficients;
                            }
                        }
                    }

                }
            );
        }

        pool.wait();

        // end second pass
        std::cout << "FOR p = " << p << " AND s = " << s << ", MAX g = " << max_g << std::endl;
        // std::cout << "EG: ";
        // for (auto coeff : max_g_v) {
        //     std::cout << coeff << ",";
        // }
        // std::cout << std::endl;

        // std::cout << "SEARCHING " << total_combinations << " COMBINATIONS OF " << quadratics.size() << " QUADRATICS" << std::endl;

        if (max_g < 3) break;

        for (size_t a = 0; a < p; a++) {
            if (singles[a] == s + 1 && quadratics.erase(a) > 0) {
                // I THINK THIS WORKS?
                std::erase(quadratic_coefficients, a);
                for (auto& subset_group : girths) {
                    std::erase_if(subset_group, [a](const std::vector<size_t>& subset) {
                        return std::binary_search(subset.begin(), subset.end(), a);
                        }
                    );
                }
                // TODO: can also erase all quadratics that will always have their girth induced <= s+1
                //.. but is it possible to actually do the removal and benefit from it?
            }
        }
    }
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

// todo: static thread pool structure initialized at start of main, used in various tasks

int main(int, char**) {

    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    // mpz_class p;
    // p = 2;
    // for (int i = 0; i < 100; i++) {
    //     compute_girth_stats(p.get_ui());
    //     mpz_nextprime(p.get_mpz_t(), p.get_mpz_t());
    //     std::cout << std::endl;
    // }

    // size_t p = 401;
    //97 = 17 seconds with twos
    //401 = 5 seconds with full lookup

    compute_girth_stats(601);



    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    std::cerr << "Finished in: " << duration << " ms" << std::endl;
}
