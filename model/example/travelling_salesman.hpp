#ifndef GENETIC_ALGORITHM_TRAVELLING_SALESMAN_HPP
#define GENETIC_ALGORITHM_TRAVELLING_SALESMAN_HPP

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <array>
#include <cmath>
#include "model/genetic_algorithm.hpp"

[[maybe_unused]] auto travelling_salesman_problem() -> void {
    using namespace ai::gal;

    auto earth_distance = [](long double lat1, long double long1, long double lat2, long double long2) -> long double {
        auto rad = [](long double deg) { return deg * ((M_PI) / 180.0);  };

        auto r_lat1 = rad(lat1);
        auto r_long1 = rad(long1);
        auto r_lat2 = rad(lat2);
        auto r_long2 = rad(long2);

        auto delta_long = r_long2 - r_long1;
        auto delta_lat = r_lat2 - r_lat1;

        auto ans = pow(sin(delta_lat / 2), 2) + cos(r_lat1) * cos(r_lat2) * pow(sin(delta_long / 2), 2);

        return 2 * asin(sqrt(ans)) * 6371l;
    };

    using Point = std::pair<double, double>;
    using ChrType = std::vector<std::size_t>;

    // Cities coordinates vector, instance from burma24.tsp
    std::vector<Point> cities {
            { 16.47, 96.10 },
            { 16.47, 94.44 },
            { 20.09, 92.54 },
            { 22.39, 93.37 },
            { 25.23, 97.24 },
            { 22.00, 96.05 },
            { 20.47, 97.02 },
            { 17.20, 96.29 },
            { 16.30, 97.38 },
            { 14.05, 98.12 },
            { 16.53, 97.38 },
            { 21.52, 95.59 },
            { 19.41, 97.13 },
            { 20.09, 94.55 }
    };

    // Lambda to calculate accumulated distance from a path
    auto accumulated_distance = [&](auto chr) {
        auto distance_at = [&](auto i, auto j) {
            auto [ lat1, long1 ] = cities.at(i);
            auto [ lat2, long2 ] = cities.at(j);
            return earth_distance(lat1, long1, lat2, long2); // or std::sqrt(std::pow(x2 - x1, 2) + std::pow(y2 - y1, 2));
        };

        return std::transform_reduce(std::begin(chr), std::prev(std::end(chr)),
                                     std::next(std::begin(chr)), 0.0l, std::plus{}, distance_at);
    };

    // Initialize random population
    auto population = random_population<ChrType>(50, [&, len = std::size(cities)](){
        return random_chromosome<typename ChrType::value_type>(len, [len]() {
            return sto::random_between(0ul, len);
        });
    });

    // Create Genetic Operators
    Best_Selector<ChrType> best_selector (1);
    Tournament_Selector<ChrType> selector (2, 0.8);
    Binary_Mating<ChrType> crossover (order_crossover<ChrType>, 0.8);
    Mutator<ChrType> mutator (swap_mutation<ChrType>, 0.02);

    // Fitness is calculated dividing one by the accumulated distance cost from the chromosome path,
    // so the longest path has the least fitness value.
    Fitness_Applier<ChrType> fitness([accumulated_distance](ChrType const& chr) {
        return 1 / accumulated_distance(chr);
    });

    // Create Genetic Evolution Context
    // Phase 1 happens first and only once per iteration.
    // * Selects all individuals and apply fitness function.

    // Phase 2 happens second and only once per iteration, often pick some organism
    // * Selects some individuals. In this case, it selects only the best one from previous iteration.

    // Phase 3 happens until the new population reaches the size of the previous one.
    // * Selects some individuals, maybe applies crossover, possibly mutation, and finally updates fitness.
    Genetic_Context gc (population, { fitness.fun() }, { best_selector.fun() },
                { Combined_Operator<ChrType>(selector, crossover, mutator, fitness).fun() });

    // Evolve population for 100 epochs, then sort it by its fitness.
    auto evolved_pop = gc.evolve_n(100).m_current_population;
    std::sort(std::begin(evolved_pop), std::end(evolved_pop), [](auto const& l, auto const& r) { return l.fitness > r.fitness; });

    for (auto const& [ chromo, fitness ] : evolved_pop) {
        for (auto c : chromo) {
            std::cout << std::setw(2) << c << ' ';
        }
        std::cout << std::setprecision(10) << fitness << ' '
                  << std::setprecision(10) << accumulated_distance(chromo) << '\n';
    }
}

#endif //GENETIC_ALGORITHM_TRAVELLING_SALESMAN_HPP
