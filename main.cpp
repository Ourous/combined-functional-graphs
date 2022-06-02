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

// TODO: look into linear-time set merging since they're all ordered

uintmax_t binomial(uintmax_t n, uintmax_t k) {
    if (n == k) return 1;
    if (n < k) return 0;

    mpz_class out;
    mpz_bin_uiui(out.get_mpz_t(), n, k);

    if (!out.fits_ulong_p()) {
        std::cerr << "Cannot fit binomial into uintmax_t: CRASHING!" << std::endl;
        exit(1);
    }

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

// other than singles, only stores quadratic properties if the single function graph has a girth > 2
struct quadratic_data {
    std::map<size_t, std::map<size_t, std::set<size_t>>> quadratics; // input-output map of quadratics in F_p
    std::vector<size_t> quadratic_coefficients; // ordered lookup of indices to coefficients in quadratics (used for combination enumeration)
    std::map<size_t, size_t> singles; // girths of single-function graphs of the corresponding coefficient
    std::map<unsigned int, std::vector<std::vector<size_t>>> girths; // girths of girth-bounding combined functional graphs
};

quadratic_data init_quad_data(size_t p, thread_pool& pool) {
    quadratic_data data;

    {
        std::mutex pass_mutex;

        const auto quadratics_per_job = p / pool.threads;
        const auto leftover_quadratics = p % pool.threads;
        size_t start_quad_id = 0;

        std::set<size_t> sorted_quadratic_coefficients;

        for (size_t which_thread = 0; which_thread < pool.threads; which_thread++) {
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
                            size_t y = (x * x + x + a) % p;
                            if (x == y) girth_is_one = true;
                            else if (!girth_maybe_two) {
                                if (x == (y * y + y + a) % p) girth_maybe_two = true;
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
                    data.singles.insert(local_singles.begin(), local_singles.end());
                    data.quadratics.insert(local_quadratics.begin(), local_quadratics.end());
                    sorted_quadratic_coefficients.insert(local_quadratic_coefficients.begin(), local_quadratic_coefficients.end());
                }
            );
            start_quad_id += quadratics_in_job;
        }

        pool.wait();

        data.quadratic_coefficients.insert(data.quadratic_coefficients.begin(), sorted_quadratic_coefficients.begin(), sorted_quadratic_coefficients.end());

    }

    {

        std::mutex pass_mutex;
        std::vector<std::pair<size_t, std::map<size_t, std::set<size_t>>>> quadratics_to_check;
        for (const auto& [a, points] : data.quadratics) {
            quadratics_to_check.emplace_back(std::pair{ a, points });
        }
        const auto quadratics_per_job = data.quadratics.size() / pool.threads;
        const auto leftover_quadratics = data.quadratics.size() % pool.threads;
        size_t start_quad_id = 0;

        for (size_t which_thread = 0; which_thread < pool.threads; which_thread++) {
            const auto quadratics_in_job = quadratics_per_job + ((leftover_quadratics > which_thread) ? 1 : 0);
            pool.do_job(
                [&, quadratics_in_job, start_quad_id] {

                    std::unordered_map<size_t, size_t> local_singles;
                    std::map<unsigned int, std::vector<std::vector<size_t>>> local_girths;

                    for (size_t quadratic_id = 0; quadratic_id < quadratics_in_job; quadratic_id++) {
                        const auto& quadratic = quadratics_to_check[quadratic_id + start_quad_id];
                        const auto a = quadratic.first;
                        const auto& points = quadratic.second;
                        csr_matrix m(p, points);
                        unsigned int g = (unsigned int)p;
                        for (size_t v = 0; v < p; v++) {
                            g = std::min(g, bfs(m, v, g));
                        }
                        local_singles[a] = g;
                        local_girths[g].emplace_back(std::vector{ a });
                    }
                    auto lock = std::scoped_lock(pass_mutex);
                    data.singles.insert(local_singles.begin(), local_singles.end());
                    for (auto& [g, subsets] : local_girths) {
                        data.girths[g].insert(data.girths[g].end(), subsets.begin(), subsets.end());
                    }
                }
            );
            start_quad_id += quadratics_in_job;
        }

        pool.wait();
    }

    return data;
}

void find_girth_bounding_cycles(size_t p) {

    std::cout << "p = " << p << std::endl;

    // initialize the various function sets
    std::map<size_t, std::map<size_t, std::set<size_t>>> quadratics;

    // quadratic_coefficients[index] = coefficient (in order)
    std::vector<size_t> quadratic_coefficients;

    // std::vector<size_t> singles(p * p * p, 0); // a -> g
    std::map<size_t, size_t> singles;
    // [girth : vector of sets of n-vectors in order of increasing n ]
    std::map<unsigned int, std::vector<std::vector<size_t>>> girths;

    size_t num_threads = std::thread::hardware_concurrency();
    thread_pool pool(num_threads);

    {
        auto data = init_quad_data(p, pool);
        quadratics = std::move(data.quadratics);
        quadratic_coefficients = std::move(data.quadratic_coefficients);
        singles = std::move(data.singles);
        girths = std::move(data.girths);
    }

    for (size_t i = 0; i < p; i++) {
        std::cout << i << ", g = " << singles[i] << std::endl;
    }

    if (girths.size() == 0) return;

    for (size_t a = 0; a < p; a++) {
        if (singles[a] < 3 && quadratics.erase(a) > 0) {
            std::erase(quadratic_coefficients, a);
            for (auto& [g, subsets] : girths) {
                std::erase_if(
                    subsets,
                    [a](const std::vector<size_t>& subset) {
                        return std::binary_search(subset.begin(), subset.end(), a);
                    }
                );

            }
        }
    }

    std::mutex girths_mutex;
    std::queue<std::pair<unsigned int, std::vector<size_t>>> coefficient_queue_2nd_pass;

    for (size_t s = 2; s <= quadratics.size() && s < p; s++) {

        std::unordered_set<size_t> coefficients_to_keep;

        // start 1st pass threads here
        unsigned int max_g = 0;
        {
            const auto total_combinations = binomial(quadratics.size(), s);
            const auto combinations_per_job = total_combinations / num_threads;
            const auto leftover_combinations = total_combinations % num_threads;

            size_t start_combination_id = 0;

            for (size_t which_thread = 0; which_thread < num_threads; which_thread++) {

                const auto combinations_in_job = combinations_per_job + ((leftover_combinations > which_thread) ? 1 : 0);

                if (combinations_in_job > 0) pool.do_job(
                    // I think the order of variable declaration and memory layout in this function is *really* important
                    //.. and should be analyzed further: moving the declaration of coefficients around seems to change the runtime by a whole second for 601
                    [&, start_combination_id, combinations_in_job] {

                        unsigned int local_max_g = 0;

                        std::vector<bool> combination_indices(quadratics.size(), false);

                        for (auto i : find_k_combination(s, start_combination_id)) combination_indices[i] = true;

                        std::vector<size_t> coefficients(s);

                        for (unsigned long combination_id = 0; combination_id < combinations_in_job; combination_id++, std::next_permutation(combination_indices.rbegin(), combination_indices.rend())) {
                            size_t coefficient_index = 0;
                            for (size_t i = 0; i < combination_indices.size(); i++) {
                                if (!combination_indices[i]) continue;
                                coefficients[coefficient_index++] = quadratic_coefficients[i];
                            }

                            // TODO: try special-casing checking for size-2 subsets, since that's most of the work being done

                            unsigned int g = 0;
                            [&] {
                                for (const auto& [g_, subsets] : girths) {
                                    // if we order the subsets lexicographically we can perform a binary search for the range in which subsets that could be inside
                                    //.. our coefficients can exist in with reference to the first and last elements of our coefficients
                                    //.. eg we can consider only the range [subset.begin() >= coefficients.begin(), subset.end() <= coefficients.end()]
                                    for (const auto& subset : subsets) {
                                        if (std::includes(coefficients.begin(), coefficients.end(), subset.begin(), subset.end())) {
                                            g = g_;
                                            return;
                                        }
                                    }
                                }
                            }();

                            if (g > s) {
                                auto lock = std::scoped_lock(girths_mutex);
                                coefficient_queue_2nd_pass.emplace(std::pair{ g, coefficients });
                            }
                            else if (g > local_max_g) local_max_g = g;

                        }

                        {
                            auto lock = std::scoped_lock(girths_mutex);
                            if (local_max_g > max_g) max_g = local_max_g;
                        }
                    }
                );

                start_combination_id += combinations_in_job;
            }


            pool.wait();
        }

        // do second pass here
        // redistribute queues

        for (size_t which_t = 0; which_t < num_threads; which_t++) {
            pool.do_job(
                [&] {
                    std::vector<size_t> coefficients;
                    unsigned int g;
                    while (true) {
                        {
                            auto lock = std::scoped_lock(girths_mutex);

                            if (coefficient_queue_2nd_pass.empty()) return;

                            auto arguments = std::move(coefficient_queue_2nd_pass.front());
                            coefficient_queue_2nd_pass.pop();
                            coefficients = std::move(arguments.second);
                            g = arguments.first;
                        }


                        std::map<size_t, std::set<size_t>> combined_function;
                        for (auto a : coefficients) {
                            for (const auto& [x, y] : quadratics[a]) {
                                if (combined_function.contains(x)) combined_function[x].insert(y.begin(), y.end());
                                else combined_function.emplace(std::pair{ x, y });
                            }
                        }

                        csr_matrix m(p, std::move(combined_function));
                        unsigned int initial_g = g;

                        for (size_t v = 0; v < p && g > s; v++) g = std::min(g, bfs(m, v, g));
                        if (initial_g > g) {
                            {
                                auto lock = std::scoped_lock(girths_mutex);

                                for (auto coeff : coefficients) {
                                    std::cout << coeff << ", ";
                                }
                                std::cout << "g = " << g << std::endl;

                                girths[g].emplace_back(std::move(coefficients));

                                if (g > max_g) {
                                    max_g = g;
                                }
                                if (g >= s) {
                                    coefficients_to_keep.insert(coefficients.begin(), coefficients.end());
                                }
                            }
                        }
                        else {
                            auto lock = std::scoped_lock(girths_mutex);
                            if (g > max_g) {
                                max_g = g;
                            }
                            if (g > s) { // this should always be true here
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

        for (size_t a = 0; a < p; a++) {

            if ((singles[a] == s + 1 || !coefficients_to_keep.contains(a)) && quadratics.erase(a) > 0) {
                std::erase(quadratic_coefficients, a);
                for (auto& [g, subsets] : girths) {
                    std::erase_if(
                        subsets,
                        [a](const std::vector<size_t>& subset) {
                            return std::binary_search(subset.begin(), subset.end(), a);
                        }
                    );

                }
            }
        }

        if (max_g <= s) break; // found all cycles
    }
}

// assumes g_limit >= 3
void compute_girth_stats(size_t p, size_t s_limit, unsigned int g_limit) {
    // initialize the various function sets
    std::map<size_t, std::map<size_t, std::set<size_t>>> quadratics;

    // quadratic_coefficients[index] = coefficient (in order)
    std::vector<size_t> quadratic_coefficients;

    // std::vector<size_t> singles(p * p * p, 0); // a -> g
    std::map<size_t, size_t> singles;
    // [girth : vector of sets of n-vectors in order of increasing n ]
    std::map<unsigned int, std::vector<std::vector<size_t>>> girths;

    size_t num_threads = std::thread::hardware_concurrency();
    thread_pool pool(num_threads);

    {
        auto data = init_quad_data(p, pool);
        quadratics = std::move(data.quadratics);
        quadratic_coefficients = std::move(data.quadratic_coefficients);
        singles = std::move(data.singles);
        girths = std::move(data.girths);
    }

    if (girths.size() == 0) return;

    for (size_t a = 0; a < p; a++) {
        if (singles[a] < g_limit && quadratics.erase(a) > 0) {
            std::erase(quadratic_coefficients, a);
            for (auto& [g, subsets] : girths) {
                std::erase_if(
                    subsets,
                    [a](const std::vector<size_t>& subset) {
                        return std::binary_search(subset.begin(), subset.end(), a);
                    }
                );

            }
        }
    }

    std::mutex girths_mutex;
    std::queue<std::pair<unsigned int, std::vector<size_t>>> coefficient_queue_2nd_pass;

    for (size_t s = 2; s <= quadratics.size() && s < s_limit; s++) {

        std::unordered_set<size_t> coefficients_to_keep;

        // start 1st pass threads here
        unsigned int max_g = 0;
        {
            const auto total_combinations = binomial(quadratics.size(), s);
            const auto combinations_per_job = total_combinations / num_threads;
            const auto leftover_combinations = total_combinations % num_threads;

            size_t start_combination_id = 0;

            for (size_t which_thread = 0; which_thread < num_threads; which_thread++) {

                const auto combinations_in_job = combinations_per_job + ((leftover_combinations > which_thread) ? 1 : 0);

                if (combinations_in_job > 0) pool.do_job(
                    // I think the order of variable declaration and memory layout in this function is *really* important
                    //.. and should be analyzed further: moving the declaration of coefficients around seems to change the runtime by a whole second for 601
                    [&, start_combination_id, combinations_in_job] {

                        std::unordered_set<size_t> local_coefficients_to_keep;

                        unsigned int local_max_g = 0;

                        std::vector<bool> combination_indices(quadratics.size(), false);

                        for (auto i : find_k_combination(s, start_combination_id)) { // TODO: an alternative method of doing this with better linear performance
                            combination_indices[i] = true;
                        }

                        std::vector<size_t> coefficients(s);

                        for (unsigned long combination_id = 0; combination_id < combinations_in_job; combination_id++, std::next_permutation(combination_indices.rbegin(), combination_indices.rend())) {
                            size_t coefficient_index = 0;
                            for (size_t i = 0; i < combination_indices.size(); i++) {
                                if (!combination_indices[i]) continue;
                                coefficients[coefficient_index++] = quadratic_coefficients[i];
                            }

                            unsigned int g = 0;
                            [&] {
                                for (const auto& [g_, subsets] : girths) {
                                    // if we order the subsets lexicographically we can perform a binary search for the range in which subsets that could be inside
                                    //.. our coefficients can exist in with reference to the first and last elements of our coefficients
                                    //.. eg we can consider only the range [subset.begin() >= coefficients.begin(), subset.end() <= coefficients.end()]
                                    for (const auto& subset : subsets) {
                                        if (std::includes(coefficients.begin(), coefficients.end(), subset.begin(), subset.end())) { // there is an idea that if I find every possible combination with a coefficient already, I can discard the coefficient
                                            g = g_;
                                            return;
                                        }
                                    }
                                }
                            }();

                            if (g <= s) {
                                if (g > local_max_g) local_max_g = g;
                                if (g >= g_limit) local_coefficients_to_keep.insert(coefficients.begin(), coefficients.end());
                            }
                            else {
                                auto lock = std::scoped_lock(girths_mutex);
                                coefficient_queue_2nd_pass.emplace(std::pair{ g, coefficients });
                            }
                        }

                        {
                            auto lock = std::scoped_lock(girths_mutex);

                            if (local_max_g > max_g) max_g = local_max_g;

                            coefficients_to_keep.insert(local_coefficients_to_keep.begin(), local_coefficients_to_keep.end());
                        }
                    }
                );

                start_combination_id += combinations_in_job;
            }


            pool.wait();
        }
        // do second pass here
        // redistribute queues

        for (size_t which_t = 0; which_t < num_threads; which_t++) {

            pool.do_job(
                [&] {
                    std::vector<size_t> coefficients;
                    unsigned int g;
                    while (true) {
                        {
                            auto lock = std::scoped_lock(girths_mutex);

                            if (coefficient_queue_2nd_pass.empty()) return;

                            auto arguments = std::move(coefficient_queue_2nd_pass.front());
                            coefficient_queue_2nd_pass.pop();
                            coefficients = std::move(arguments.second);
                            g = arguments.first;
                        }

                        std::map<size_t, std::set<size_t>> combined_function;
                        for (auto a : coefficients) {
                            for (const auto& [x, y] : quadratics[a]) {
                                if (combined_function.contains(x)) combined_function[x].insert(y.begin(), y.end());
                                else combined_function.emplace(std::pair{ x, y });
                            }
                        }

                        csr_matrix m(p, std::move(combined_function));
                        unsigned int initial_g = g;

                        for (size_t v = 0; v < p && g > s; v++) g = std::min(g, bfs(m, v, g));
                        if (initial_g > g) { // NOTE: there's some strange performance impact depending on whether I lock before or after the branch
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
                            }
                            if (g >= g_limit) {
                                coefficients_to_keep.insert(coefficients.begin(), coefficients.end());
                            }
                        }
                    }

                }
            );
        }

        pool.wait();

        // end second pass

        if (max_g < g_limit) {
            std::cout << "p = " << p << ", s = " << s << ", max g < " << g_limit << std::endl;
            return;
        };

        std::cout << "p = " << p << ", s = " << s << ", max g = " << max_g << std::endl;

        for (size_t a = 0; a < p; a++) {

            if (!coefficients_to_keep.contains(a) && quadratics.erase(a) > 0) {

                std::erase(quadratic_coefficients, a);
                for (auto& [g, subsets] : girths) {
                    std::erase_if(
                        subsets,
                        [a](const std::vector<size_t>& subset) {
                            return std::binary_search(subset.begin(), subset.end(), a);
                        }
                    );

                }
            }
        }

        if (quadratics.size() < s + 1) {
            std::cout << "p = " << p << ", s = " << s + 1 << ", max g < " << g_limit << std::endl;
            return;
        }
    }
}

int main(int, char**) {

    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    mpz_class p;
    p = 0;
    while (p.get_ui() < 3000) {
        mpz_nextprime(p.get_mpz_t(), p.get_mpz_t());
        // compute_girth_stats(p.get_ui(), p.get_ui(), 3);
        find_girth_bounding_cycles(p.get_ui());
        std::cout << std::endl;
    }

    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    std::cerr << "Finished in: " << duration << " ms" << std::endl;
}
