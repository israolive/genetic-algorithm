#ifndef GENETIC_ALGORITHM_NQUEENS_HPP
#define GENETIC_ALGORITHM_NQUEENS_HPP

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <array>
#include <cmath>
#include "model/genetic_algorithm.hpp"

[[maybe_unused]] auto nqueens_problem() -> void {
    using namespace ai::gal;

    using Point = struct { long x, y; };

    using ChrType = std::vector<long>;
    using AlleleType = typename ChrType::value_type;

    auto n = 8ull;

    auto can_kill_occur = [](auto const& p1, auto const& p2) {
        return p1.x == p2.x || p1.y == p2.y || std::abs(p2.x - p1.x) == std::abs(p2.y - p1.y);
    };

    auto count_kills = [=](auto const &queens) {
        auto alive = std::vector<bool>(n, true);

        std::for_each(std::begin(queens), std::end(queens), [&, queen_i = 0l](auto const killer) mutable {
            std::transform(std::begin(queens), std::end(queens), std::begin(alive), std::begin(alive),
                           [&, queen_j = 0l](auto const target, auto const lives) mutable {
                               return lives && !(queen_i != queen_j &&
                                        can_kill_occur(Point{queen_i, killer}, Point{queen_j++, target}));
                           });
            ++queen_i;
        });

        return std::count(std::begin(alive), std::end(alive), false);
    };

    auto chr_kill_factor = [count_kills, n](auto const& chr) {
        return 1 - (count_kills(chr) / static_cast<double>(n));
    };

    // Initialize random population
    auto population = random_population<ChrType>(50, [n](){
        return random_chromosome<AlleleType>(n, [n]() {
            return sto::random_between(0l, static_cast<AlleleType>(n));
        });
    });

    // Create Genetic Operators
    auto best_selector = Best_Selector<ChrType>(1);
    auto selector = Tournament_Selector<ChrType>(2, 0.8);
    auto crossover = Binary_Mating<ChrType>(order_crossover<ChrType>, 0.8);
    auto mutator = Mutator<ChrType>(swap_mutation<ChrType>, 0.02);
    auto fitness = Fitness_Applier<ChrType>(chr_kill_factor);

    Genetic_Context gc (population,{ fitness.fun() }, { best_selector.fun() },
                        { Combined_Operator<ChrType>(selector, crossover, mutator, fitness).fun() });

    // Evolve population for 100 epochs, then sort it by its fitness.
    auto evolved_pop = gc.evolve_n(1000).m_current_population;
    std::sort(std::begin(evolved_pop), std::end(evolved_pop), [](auto const& l, auto const& r) { return l.fitness > r.fitness; });

    for (auto const& [ chromo, fitness ] : evolved_pop) {
        for (auto c : chromo) {
            std::cout << std::setw(2) << c << ' ';
        }
        std::cout << std::setprecision(10) << fitness << ' '
                  << std::setprecision(10) << count_kills(chromo) << '\n';
    };
}

#endif //GENETIC_ALGORITHM_NQUEENS_HPP
