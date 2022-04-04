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
#include <unordered_map>
#include <unordered_set>

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
    if (!out.fits_ulong_p()) {
        std::cout << "Cannot fit binomial into uintmax_t: CRASHING!" << std::endl;
        exit(1);
    }
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

void compute_girth_stats(size_t p, size_t s_limit, unsigned int g_limit) {
    // initialize the various function sets
    std::map<size_t, std::map<size_t, std::set<size_t>>> quadratics;

    // quadratic_coefficients[index] = coefficient (in order)
    std::vector<size_t> quadratic_coefficients;

    // std::vector<size_t> singles(p * p * p, 0); // a -> g
    std::unordered_map<size_t, size_t> singles;
    // [girth : vector of sets of n-vectors in order of increasing n ]
    std::map<unsigned int, std::vector<std::vector<size_t>>> girths;

    // TODO: add multithreading to this first part

    size_t t = std::thread::hardware_concurrency();
    thread_pool pool(t);

    // for (size_t a = 0; a < p; a++) {
    //     std::map<size_t, std::set<size_t>> points;
    //     bool girth_is_one = false, girth_maybe_two = false;
    //     for (size_t x = 0; x < p && !girth_is_one; x++) { // works for p < 2^63
    //         size_t y = (x * x + a) % p;
    //         if (x == y) girth_is_one = true;
    //         else if (!girth_maybe_two) {
    //             if (x == (y * y + a) % p) girth_maybe_two = true;
    //             else points[x].insert(y);
    //         }
    //     }
    //     if (girth_is_one) {
    //         singles[a] = 1;
    //     }
    //     else if (girth_maybe_two) {
    //         singles[a] = 2;
    //     }
    //     else {
    //         quadratics.insert({ a, points });
    //         quadratic_coefficients.push_back(a);
    //     }
    // }

    {
        std::mutex pass_mutex;

        const auto quadratics_per_job = p / t;
        const auto leftover_quadratics = p % t;
        size_t start_quad_id = 0;

        std::set<size_t> sorted_quadratic_coefficients;

        for (size_t which_thread = 0; which_thread < t; which_thread++) {
            const auto quadratics_in_job = quadratics_per_job + ((leftover_quadratics > which_thread) ? 1 : 0);
            pool.do_job(
                [&, quadratics_in_job, start_quad_id] {

                    std::unordered_map<size_t, size_t> local_singles;
                    std::vector<size_t> local_quadratic_coefficients;
                    std::map<size_t, std::map<size_t, std::set<size_t>>> local_quadratics;

                    for (size_t quadratic_id = 0; quadratic_id < quadratics_in_job; quadratic_id++) {
                        auto a = quadratic_id + start_quad_id;
                        std::map<size_t, std::set<size_t>> points;
                        bool girth_is_one = false, girth_maybe_two = false;
                        for (size_t x = 0; x < p && !girth_is_one; x++) { // works for p < 2^63
                            size_t y = (x * x + a) % p;
                            if (x == y) girth_is_one = true;
                            else if (!girth_maybe_two) {
                                if (x == (y * y + a) % p) girth_maybe_two = true;
                                else points[x].insert(y);
                            }
                        }
                        if (girth_is_one) {
                            local_singles[a] = 1;
                        }
                        else if (girth_maybe_two) {
                            local_singles[a] = 2;
                        }
                        else {
                            local_quadratics.insert({ a, points });
                            local_quadratic_coefficients.push_back(a);
                        }
                    }
                    auto lock = std::scoped_lock(pass_mutex);
                    singles.insert(local_singles.begin(), local_singles.end());
                    quadratics.insert(local_quadratics.begin(), local_quadratics.end());
                    sorted_quadratic_coefficients.insert(local_quadratic_coefficients.begin(), local_quadratic_coefficients.end());
                }
            );
            start_quad_id += quadratics_in_job;
        }

        pool.wait();

        quadratic_coefficients.insert(quadratic_coefficients.begin(), sorted_quadratic_coefficients.begin(), sorted_quadratic_coefficients.end());

    }

    // for (size_t c = 1; c < p; c++) for (size_t b = 0; b < p; b++) for (size_t a = 0; a < p; a++) {
    //     std::map<size_t, std::set<size_t>> points;
    //     bool girth_is_one = false, girth_maybe_two = false;
    //     for (size_t x = 0; x < p && !girth_is_one; x++) { // works for p < 2^63
    //         size_t y = (c * x * x + b * x + a) % p;
    //         if (x == y) girth_is_one = true;
    //         else if (!girth_maybe_two) {
    //             if (x == (c * y * y + b * y + a) % p) girth_maybe_two = true;
    //             else points[x].insert(y);
    //         }
    //     }
    //     if (girth_is_one) {
    //         singles[c * p * p + b * p + a] = 1;
    //     }
    //     else if (girth_maybe_two) {
    //         singles[c * p * p + b * p + a] = 2;
    //     }
    //     else {
    //         quadratics.insert({ c * p * p + b * p + a, points });
    //         quadratic_coefficients.push_back(c * p * p + b * p + a);
    //     }
    // }

    {

        std::mutex pass_mutex;
        std::vector<std::pair<size_t, std::map<size_t, std::set<size_t>>>> quadratics_to_check;
        for (const auto& [a, points] : quadratics) {
            quadratics_to_check.emplace_back(std::pair{ a, points });
        }
        const auto quadratics_per_job = quadratics.size() / t;
        const auto leftover_quadratics = quadratics.size() % t;
        size_t start_quad_id = 0;
        // auto start_quad_iter = quadratics.begin();
        for (size_t which_thread = 0; which_thread < t; which_thread++) {
            const auto quadratics_in_job = quadratics_per_job + ((leftover_quadratics > which_thread) ? 1 : 0);
            pool.do_job(
                [&, quadratics_in_job, start_quad_id] {

                    std::unordered_map<size_t, size_t> local_singles;
                    std::map<unsigned int, std::vector<std::vector<size_t>>> local_girths;

                    for (size_t quadratic_id = 0; quadratic_id < quadratics_in_job; quadratic_id++) {
                        auto quadratic = quadratics_to_check[quadratic_id + start_quad_id];
                        auto a = quadratic.first;
                        auto points = quadratic.second;
                        csr_matrix m(p, points);
                        unsigned int g = (unsigned int)p;
                        for (size_t v = 0; v < p; v++) {
                            g = std::min(g, bfs(m, v, g));
                        }
                        local_singles[a] = g;
                        local_girths[g].emplace_back(std::vector{ a });
                    }
                    auto lock = std::scoped_lock(pass_mutex);
                    singles.insert(local_singles.begin(), local_singles.end());
                    girths.insert(local_girths.begin(), local_girths.end());
                }
            );
            start_quad_id += quadratics_in_job;
        }

        pool.wait();
    }

    // for (auto [a, points] : quadratics) {
    //     csr_matrix m(p, points);
    //     unsigned int g = (unsigned int)p;
    //     for (size_t v = 0; v < p; v++) {
    //         g = std::min(g, bfs(m, v, g));
    //     }
    //     singles[a] = g;
    //     girths[g].emplace_back(std::vector{ a });
    //     // std::cout << a << ": " << g << std::endl;
    // }
    if (girths.size() == 0) return;
    std::cout << "p = " << p << ", s = 1, max g = " << girths.rbegin()->first << std::endl;

    for (size_t a = 0; a < p; a++) {
        if (singles[a] < g_limit && quadratics.erase(a) > 0) {
            // I THINK THIS WORKS?
            std::erase(quadratic_coefficients, a);
            for (auto& [g, subsets] : girths) {
                std::erase_if(
                    subsets,
                    [a](const std::vector<size_t>& subset) {
                        return std::binary_search(subset.begin(), subset.end(), a);
                    }
                );

            }
            // TODO: can also erase all quadratics that will always have their girth induced by a subset
            //.. but is it possible to actually do the removal and benefit from it?
            // ANSWER: yes, try this next and see if we get any improvement
        }
    }



    std::mutex girths_mutex;
    std::queue<std::pair<unsigned int, std::vector<size_t>>> coefficient_queue_2nd_pass;
    size_t s = 2;
    for (; s <= quadratics.size() && s < s_limit; s++) {

        std::unordered_set<size_t> coefficients_to_keep;

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

        const auto total_combinations = binomial(quadratics.size(), s);
        const auto combinations_per_job = total_combinations / t;
        const auto leftover_combinations = total_combinations % t;

        size_t start_combination_id = 0;

        for (size_t which_thread = 0; which_thread < t; which_thread++) {

            const auto combinations_in_job = combinations_per_job + ((leftover_combinations > which_thread) ? 1 : 0);

            if (combinations_in_job > 0) pool.do_job(
                // I think the order of variable declaration and memory layout in this function is *really* important
                //.. and should be analyzed further: moving the declaration of coefficients around seems to change the runtime by a whole second for 601
                [&, start_combination_id, combinations_in_job] {

                    std::unordered_set<size_t> local_coefficients_to_keep;

                    unsigned int local_max_g = 0;
                    std::vector<size_t> local_max_g_v;

                    std::vector<bool> combination_indices(quadratics.size(), false);

                    for (auto i : find_k_combination(s, start_combination_id)) { // TODO: just use k-combination index enumeration instead, it's faster
                        combination_indices[i] = true;
                    }

                    std::vector<size_t> coefficients(s);

                    for (unsigned long combination_id = 0; combination_id < combinations_in_job; combination_id++, std::next_permutation(combination_indices.rbegin(), combination_indices.rend())) {
                        size_t coefficient_index = 0;
                        for (size_t i = 0; i < combination_indices.size(); i++) {
                            if (!combination_indices[i]) continue;
                            // coefficients.push_back(quadratic_coefficients[i]);
                            coefficients[coefficient_index++] = quadratic_coefficients[i];
                        }
                        // memory locality of these subsets may have a massive performance impact

                        // TODO: try special-casing checking for size-2 subsets with a fancy algorithm, since that's most of the work being done

                        unsigned int g = 0;
                        [&] {
                            for (const auto& [g_, subsets] : girths) {
                                // if we order the subsets lexicographically we can perform a binary search for the range in which subsets that could be inside
                                //.. our coefficients can exist in with reference to the first and last elements of our coefficients
                                //.. eg we can consider only the range [subset.begin() >= coefficients.begin(), subset.end() <= coefficients.end()]
                                for (const auto& subset : subsets) {
                                    // TODO: check if std::includes starts with a binary search for the elements at the ends of the ranges
                                    // maybe if I sort `subset` lexicographically?
                                    if (std::includes(coefficients.begin(), coefficients.end(), subset.begin(), subset.end())) { // there is an idea that if I find every possible combination with a coefficient already, I can discard the coefficient
                                        g = g_;
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
                            if (g >= g_limit) {
                                local_coefficients_to_keep.insert(coefficients.begin(), coefficients.end());
                            }
                        }
                        else {
                            auto lock = std::scoped_lock(girths_mutex);
                            coefficient_queue_2nd_pass.emplace(std::pair{ g, coefficients });
                            // maybe also store g here
                        }
                    }

                    {
                        auto lock = std::scoped_lock(girths_mutex);
                        if (local_max_g > max_g) {
                            max_g = local_max_g;
                            max_g_v = local_max_g_v;
                        }
                        coefficients_to_keep.insert(local_coefficients_to_keep.begin(), local_coefficients_to_keep.end());
                    }
                }
            );
            start_combination_id += combinations_in_job;
        }


        pool.wait();

        // for (auto [coeff, g] : max_girth_per_coefficient) {
        //     std::cout << coeff << ": " << g << std::endl;
        // }

        // do second pass here
        // redistribute queues

        for (size_t which_t = 0; which_t < t; which_t++) {
            // if (coefficient_queues_2nd_pass[which_t].empty()) continue;
            pool.do_job(
                [&] {
                    std::vector<size_t> coefficients;
                    unsigned int g;
                    while (true) {
                        {
                            auto lock = std::scoped_lock(girths_mutex);
                            // maybe use cond_variable here, or scoped_lock
                            if (coefficient_queue_2nd_pass.empty()) return;

                            auto arguments = std::move(coefficient_queue_2nd_pass.front());
                            coefficient_queue_2nd_pass.pop();
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
                                girths[g].emplace_back(std::move(coefficients));
                                if (g > max_g) {
                                    max_g = g;
                                }
                                if (g >= g_limit) {
                                    coefficients_to_keep.insert(coefficients.begin(), coefficients.end());
                                }
                            }
                        }
                        else {
                            auto lock = std::scoped_lock(girths_mutex);
                            if (g > max_g) {
                                max_g = g;
                                max_g_v = coefficients;
                            }
                            if (g >= g_limit) {
                                coefficients_to_keep.insert(coefficients.begin(), coefficients.end());
                            }
                        }
                    }

                }
            );
        }

        // TODO: after there are no more girth-inducing subsets to be found (max_g = s), restructure `girths` so that it's more efficient for lookups
        //.. for example, maybe I could determine which minimum and maximum coefficients a given thread would be looking at
        //.. and extract these from the main girths mapping into a flat girths vector, ordered by g~>s
        //.. and maybe even further ordered by lexicographic order (but I may have to again subdivide the data structure which might not help)

        pool.wait();

        // for (auto c : max_g_v) {
        //     std::cout << c / (p * p) << "x^2 + " << (c / p) % p << "x + " << c % p << std::endl;
        // }

        // end second pass

        if (max_g < g_limit) {
            std::cout << "p = " << p << ", s = " << s << ", max g < " << g_limit << std::endl;
            return;
        };

        std::cout << "p = " << p << ", s = " << s << ", max g = " << max_g << std::endl;

        for (size_t a = 0; a < p; a++) {

            if (!coefficients_to_keep.contains(a) && quadratics.erase(a) > 0) { // or g < 3?
                // I THINK THIS WORKS?
                std::erase(quadratic_coefficients, a);
                for (auto& [g, subsets] : girths) {
                    std::erase_if(
                        subsets,
                        [a](const std::vector<size_t>& subset) {
                            return std::binary_search(subset.begin(), subset.end(), a);
                        }
                    );

                }
                // TODO: can also erase all quadratics that will always have their girth induced by a subset
                //.. but is it possible to actually do the removal and benefit from it?
                // ANSWER: yes, try this next and see if we get any improvement
            }
        }

        // std::cout << "FOR p = " << p << " AND s = " << s << ", MAX g = " << max_g << std::endl;

        if (quadratics.size() < s + 1) {
            std::cout << "p = " << p << ", s = " << s + 1 << ", max g < " << g_limit << std::endl;
            return;
        }

        // std::cout << "FOR p = " << p << " AND s = " << s << ", MAX g = " << max_g << std::endl;

        // for (size_t a = 0; a < p; a++) {
        //     if (singles[a] == s + 1 && quadratics.erase(a) > 0) {
        //         // I THINK THIS WORKS?
        //         std::erase(quadratic_coefficients, a);
        //         for (auto& [g, subsets] : girths) {
        //             std::erase_if(
        //                 subsets,
        //                 [a](const std::vector<size_t>& subset) {
        //                     return std::binary_search(subset.begin(), subset.end(), a);
        //                 }
        //             );

        //         }
        //         // TODO: can also erase all quadratics that will always have their girth induced by a subset
        //         //.. but is it possible to actually do the removal and benefit from it?
        //         // ANSWER: yes, try this next and see if we get any improvement
        //     }
        // }
    //     auto initial = quadratics.size();
    // start_precog:
    //     for (size_t i = 0; i < quadratic_coefficients.size(); i++) {
    //         bool will_erase = true;
    //         unsigned int most_g = 0;
    //         for (size_t j = 0; will_erase && j < quadratic_coefficients.size(); j++) {
    //             size_t a, b;
    //             // auto k = j <=> i; // TODO: maybe use 3-way comparison
    //             if (j == i) continue;
    //             if (i < j) {
    //                 a = quadratic_coefficients[i];
    //                 b = quadratic_coefficients[j];
    //             }
    //             if (j < i) {
    //                 a = quadratic_coefficients[j];
    //                 b = quadratic_coefficients[i];
    //             }
    //             for (auto& [g, subsets] : girths) {
    //                 // if (g > s) break;
    //                 for (auto& subset : subsets) {
    //                     if (subset.size() == 1) continue;
    //                     if (subset.size() > 2) break;
    //                     if (subset[0] == a && subset[1] == b) {
    //                         most_g = std::max(g, most_g);
    //                         goto found_subset;
    //                     }
    //                 }
    //             }
    //             will_erase = false;
    //         found_subset:;
    //         }
    //         if (will_erase && most_g < 3) {
    //             size_t a = quadratic_coefficients[i];
    //             // std::cout << "for a = " << a << ", the highest g = " << most_g << std::endl;
    //             std::erase(quadratic_coefficients, a);
    //             quadratics.erase(a);
    //             for (auto& [g, subsets] : girths) {
    //                 std::erase_if(
    //                     subsets,
    //                     [a](const std::vector<size_t>& subset) {
    //                         return std::binary_search(subset.begin(), subset.end(), a);
    //                     }
    //                 );
    //             }
    //             goto start_precog;
    //         } // TODO: potentially have to re-run the previous erasure step here to get rid of subsets that don't have these quadratics in them anymore
    //     }
    //     if (quadratics.size() != initial) {
    //         std::cout << "saved " << initial - quadratics.size() << std::endl;
    //     }
       // if (max_g <= s) break; // found all cycles
    }
    /*
        std::vector<std::pair<unsigned int, std::vector<size_t>>> flat_girths;

        for (auto [g, subsets] : girths) {
            for (auto subset : subsets) {
                flat_girths.emplace_back(std::pair{ g, subset });
            }
        }

        for (; s <= quadratics.size() && s < s_limit; s++) {
            std::unordered_set<size_t> coefficients_to_keep;

            // todo: split cofficients into (t) queues that can be taken out of by the (t) threads
            //.. and either discarded by the threads or stored in another queue
            //.. when all coefficients are initially examined against the girth subset lookup
            //.. (possibly redistribute the members of these queues evenly first)
            //.. stop these (t) threads and start another (t) threads on the (t) new queues
            //.. which will perform the girth computation and store it until they're done all the work in their queue
            //.. after which they all write their results into the shared girths

            // start 1st pass threads here
            unsigned int max_g = 0;

            const auto total_combinations = binomial(quadratics.size(), s);
            const auto combinations_per_job = total_combinations / t;
            const auto leftover_combinations = total_combinations % t;

            size_t start_combination_id = 0;

            for (size_t which_t = 0; which_t < t; which_t++) {

                const auto combinations_in_job = combinations_per_job + ((leftover_combinations > which_t) ? 1 : 0);

                if (combinations_in_job > 0) pool.do_job(
                    [&, start_combination_id, combinations_in_job] {

                        std::vector<std::pair<unsigned int, std::vector<size_t>>> local_flat_girths;

                        for (auto& entry : flat_girths) {
                            if (entry.second.back() <= quadratic_coefficients[find_k_combination(s, start_combination_id + combinations_in_job - 1).front()]) {
                                local_flat_girths.push_back(entry);
                            }
                        }

                        std::unordered_set<size_t> local_coefficients_to_keep;

                        unsigned int local_max_g = 0;

                        std::vector<bool> combination_indices(quadratics.size(), false);

                        for (auto i : find_k_combination(s, start_combination_id)) { // TODO: just use k-combination index enumeration instead, it's faster
                            combination_indices[i] = true;
                        }

                        std::vector<size_t> coefficients(s);

                        for (unsigned long combination_id = 0; combination_id < combinations_in_job; combination_id++, std::next_permutation(combination_indices.rbegin(), combination_indices.rend())) {
                            size_t coefficient_index = 0;
                            for (size_t i = 0; i < combination_indices.size(); i++) {
                                if (!combination_indices[i]) continue;
                                // coefficients.push_back(quadratic_coefficients[i]);
                                coefficients[coefficient_index++] = quadratic_coefficients[i];
                            }
                            // memory locality of these subsets may have a massive performance impact

                            // TODO: try special-casing checking for size-2 subsets with a fancy algorithm, since that's most of the work being done

                            unsigned int g = 0;
                            [&] {
                                for (const auto& entry : local_flat_girths) {
                                    if (std::includes(coefficients.begin(), coefficients.end(), entry.second.begin(), entry.second.end())) { // there is an idea that if I find every possible combination with a coefficient already, I can discard the coefficient
                                        g = entry.first;
                                        return;
                                    }
                                }
                            }();

                            if (g > local_max_g) {
                                local_max_g = g;
                            }
                            if (g >= g_limit) {
                                local_coefficients_to_keep.insert(coefficients.begin(), coefficients.end());
                            }
                        }

                        {
                            auto lock = std::scoped_lock(girths_mutex);
                            if (local_max_g > max_g) {
                                max_g = local_max_g;
                            }
                            coefficients_to_keep.insert(local_coefficients_to_keep.begin(), local_coefficients_to_keep.end());
                        }
                    }
                );
                start_combination_id += combinations_in_job;
            }

            pool.wait();

            if (max_g < g_limit) {
                std::cout << "FOR p = " << p << " AND s = " << s << ", MAX g < " << g_limit << std::endl;
                return;
            };

            std::cout << "FOR p = " << p << " AND s = " << s << ", MAX g = " << max_g << std::endl;

            for (size_t a = 0; a < p; a++) {

                if (!coefficients_to_keep.contains(a) && quadratics.erase(a) > 0) { // or g < 3?
                    // I THINK THIS WORKS?
                    std::erase(quadratic_coefficients, a);
                    // for (auto& [g, subsets] : girths) {
                    //     std::erase_if(
                    //         subsets,
                    //         [a](const std::vector<size_t>& subset) {
                    //             return std::binary_search(subset.begin(), subset.end(), a);
                    //         }
                    //     );
                    // }
                    std::erase_if(
                        flat_girths,
                        [a](const std::pair<unsigned int, std::vector<size_t>>& entry) {
                            return std::binary_search(entry.second.begin(), entry.second.end(), a);
                        }
                    );
                    // TODO: can also erase all quadratics that will always have their girth induced by a subset
                    //.. but is it possible to actually do the removal and benefit from it?
                    // ANSWER: yes, try this next and see if we get any improvement
                }
            }

            // std::cout << "FOR p = " << p << " AND s = " << s << ", MAX g = " << max_g << std::endl;

            if (quadratics.size() < s + 1) {
                std::cout << "FOR p = " << p << " AND s = " << s + 1 << ", MAX g < " << g_limit << std::endl;
                return;
            }
        }*/
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
    mpz_class p;
    p = 0;
    // p = 0;
    // for (int i = 0; i < 100; i++) {
    while (p.get_ui() < 10000) {
        mpz_nextprime(p.get_mpz_t(), p.get_mpz_t());
        compute_girth_stats(p.get_ui(), p.get_ui(), 5);
        std::cout << std::endl;
    }

    // size_t p = 401;
    //97 = 17 seconds with twos
    //401 = 5 seconds with full lookup

    // compute_girth_stats(1009);

    // TODO: Also make graphs for s inducing g < 4, g < 5, etc.. to show the same pattern for larger g and larger p

    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    std::cerr << "Finished in: " << duration << " ms" << std::endl;
}
