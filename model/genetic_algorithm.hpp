#ifndef GENETICALGORITHM_HPP
#define GENETICALGORITHM_HPP

#include <algorithm>
#include <numeric>
#include <functional>
#include <iostream>
#include <vector>

#include "stochastic.hpp"

namespace ai::gal {

    template <typename Allele>
    using Chromosome = std::vector<Allele>;

    template <typename ChromoType>
    struct Organism {
        ChromoType c;
        double fitness {};
    };

    template <typename ChromoType>
    using Population = std::vector<Organism<ChromoType>>;

    template <typename ChromoType>
    using Genetic_Operator_Function = std::function<auto (Population<ChromoType> const&) -> Population<ChromoType>>;

    template <typename ChromoType>
    struct Genetic_Operator {
        virtual auto operate(Population<ChromoType> const& pop) const -> Population<ChromoType> = 0;

        // /todo maybe ub if derived class gets temp instanced then fun is called
        virtual auto fun() const -> Genetic_Operator_Function<ChromoType> {
            return [this](Population<ChromoType> const& pop) {
                return this->operate(pop);
            };
        }
    };

    template <typename ChromoType>
    auto operator|(Population<ChromoType> const& pop, Genetic_Operator<ChromoType> const& op) -> Population<ChromoType> {
        return op.operate(pop);
    }

    template <typename ChromoType>
    auto operator+(Population<ChromoType> const& l, Population<ChromoType> const& r) -> Population<ChromoType> {
        Population<ChromoType> pool;
        pool.reserve(std::size(l) + std::size(r));

        std::copy(std::begin(l), std::end(l), std::back_inserter(pool));
        std::copy(std::begin(r), std::end(r), std::back_inserter(pool));

        return pool;
    }

    template <typename AlleleType, typename AlleleGen>
    auto random_chromosome(std::size_t chr_size, AlleleGen gen) -> Chromosome<AlleleType> {
        return sto::n_unique_random<AlleleType>(chr_size, gen);
    }

    template <typename ChromoType, typename ChromoGen>
    auto random_population(std::size_t pop_size, ChromoGen gen) -> Population<ChromoType> {
        Population<ChromoType> pool;
        pool.reserve(pop_size);

        std::generate_n(std::back_inserter(pool), pop_size, [gen]() { return Organism<ChromoType>{ gen() }; });

        return pool;
    }

    template <typename ChromoType, typename ... GenOps>
    auto combine_genetic_operators(GenOps&& ... gen_ops) -> Genetic_Operator_Function<ChromoType> {
        return [&](Population<ChromoType> const& pop) {
            return (pop | ... | gen_ops);
        };
    }

    namespace {
        template <typename ChromoType>
        inline auto random_binary_tournament(Population<ChromoType> const& population, double best_chance = 0.8) -> std::size_t {
            std::size_t fst = sto::random_between(0ul, std::size(population)),
                        snd = sto::random_between(0ul, std::size(population));

            if (sto::rng() < best_chance) {
                if (population.at(fst).fitness > population.at(snd).fitness) { return fst; } else { return snd; }
            } else {
                if (population.at(fst).fitness < population.at(snd).fitness) { return fst; } else { return snd; }
            }
        }
    }

    template <typename ChromoType>
    struct Tournament_Selector : public Genetic_Operator<ChromoType> {
        std::size_t m_count;
        double m_best_chance;

        Tournament_Selector(std::size_t count, double best_chance)
            : m_count(count), m_best_chance(best_chance) {}

        auto operate(Population<ChromoType> const& pop) const -> Population<ChromoType> override {
            auto idxes = sto::n_unique_random<std::size_t>(m_count, [pop, chance = m_best_chance]() -> std::size_t {
                return random_binary_tournament(pop, chance);
            });

            Population<ChromoType> pool;
            pool.reserve(m_count);

            std::transform(std::begin(idxes), std::end(idxes), std::back_inserter(pool),
                           [&pop](auto idx) {
                               return pop.at(idx);
                           });

            return pool;
        }
    };

    template <typename ChromoType>
    auto make_tournament_selector(std::size_t count, double best_chance) -> Genetic_Operator<ChromoType> {
        Tournament_Selector<ChromoType> selector(count, best_chance);
        return selector;
    }

    template <typename ChromoType>
    struct Best_Selector : public Genetic_Operator<ChromoType> {
        std::size_t m_count;

        Best_Selector(std::size_t count)
            : m_count(count) {}

        auto operate(Population<ChromoType> const& pop) const -> Population<ChromoType> override {
            return {
                *std::max_element(std::begin(pop), std::end(pop), [](auto const& l, auto const& r) {
                    return l.fitness < r.fitness;
                })
            };
        }
    };

    template <typename ChromoType>
    struct Binary_Mating : public Genetic_Operator<ChromoType> {
        using Binary_Crossover_Fn = std::function<auto (ChromoType const&, ChromoType const&) -> std::pair<ChromoType, ChromoType>>;

        Binary_Crossover_Fn m_crossover;
        double m_probability;

        Binary_Mating(Binary_Crossover_Fn crossover_fn, double probability)
            : m_crossover(crossover_fn), m_probability(probability) {}

        auto operate(Population<ChromoType> const& pop) const -> Population<ChromoType> override {
            auto remainder = pop.size() % 2;

            Population<ChromoType> pool;
            pool.reserve(pop.size());

            // \todo replace this with an actual algorithm (zip tail), std::transform(begin, prev(end), next(begin), bin_op)
            for (std::size_t idx = 0; idx < pop.size() - remainder; idx += 2) {
                auto [ first, second ] = [x = m_crossover, p = m_probability, c1 = pop.at(idx).c, c2 = pop.at(idx + 1).c]() {
                    if (sto::rng() < p) { return x(c1, c2); }
                    else { return std::make_pair(c1, c2); }
                }();

                pool.push_back({ first });
                pool.push_back({ second });
            }

            if (remainder) {
                pool.push_back(pop.back());
            }

            return pool;
        }
    };

    template <typename ChromoType>
    struct Mutator : public Genetic_Operator<ChromoType> {
        using Mutator_Fn = std::function<auto (ChromoType const&) -> ChromoType>;

        Mutator_Fn m_mutation;
        double m_probability;

        Mutator(Mutator_Fn mutator_fn, double probability)
            : m_mutation(mutator_fn), m_probability(probability) {}

        auto operate(Population<ChromoType> const& pop) const -> Population<ChromoType> override {
            Population<ChromoType> pool;
            pool.reserve(pop.size());

            std::transform(std::begin(pop), std::end(pop), std::back_inserter(pool),
                           [mutation = m_mutation, p = m_probability](auto const& organism) {
                               if (sto::rng() < p) {
                                   return Organism<ChromoType> { mutation(organism.c) };
                               } else {
                                   return organism;
                               }
                           });

            return pool;
        }
    };

    template <typename ChromoType>
    struct Fitness_Applier : public Genetic_Operator<ChromoType> {
        using Fitness_Fn = std::function<auto (ChromoType const&) -> double>;

        Fitness_Fn m_measurement;

        explicit Fitness_Applier(Fitness_Fn measurement)
            : m_measurement(measurement) {}

        auto operator()(Fitness_Fn measurement) -> Fitness_Applier {
            Fitness_Applier applier = *this;
            applier.m_measurement = measurement;

            return applier;
        }

        auto operate(Population<ChromoType> const& pop) const -> Population<ChromoType> override {
            Population<ChromoType> pool;
            pool.reserve(std::size(pop));

            std::transform(std::begin(pop), std::end(pop), std::back_inserter(pool),
                           [measurer = m_measurement](auto const& organism) -> Organism<ChromoType> {
                               return { organism.c, measurer(organism.c) };
                           });

            return pool;
        }
    };

    template <typename ChromoType>
    struct Combined_Operator : Genetic_Operator<ChromoType> {
        Genetic_Operator_Function<ChromoType> m_fun;

        template <typename ... GenOps>
        explicit Combined_Operator(GenOps && ... gen_ops)
            : m_fun(combine_genetic_operators<ChromoType>(std::forward<GenOps>(gen_ops)...)) {}

        auto operate(Population<ChromoType> const& pop) const -> Population<ChromoType> override {
            return m_fun(pop);
        }

        auto fun() const -> Genetic_Operator_Function<ChromoType> {
            return m_fun;
        }
    };

    template <typename ChromoType>
    auto order_crossover(ChromoType const& l, ChromoType const& r) -> std::pair<ChromoType, ChromoType> {
        std::size_t rnd1 = sto::random_between(0ul, std::size(l)),
                    rnd2 = sto::random_between(0ul, std::size(l));

        while (rnd2 == rnd1) {
            rnd2 = sto::random_between(0ul, std::size(l));
        }

        auto offspring = [low = std::min(rnd1, rnd2), high = std::max(rnd1, rnd2) + 1] (auto const& p1, auto const& p2) {
            auto sibling = p1, begin_lookup = std::begin(p2);

            auto lookup_unique = [&]() mutable {
                begin_lookup = std::find_if(begin_lookup, std::end(p2), [&sibling, low, high](auto other_gene) {
                    return std::none_of(std::begin(sibling) + low, std::begin(sibling) + high, [other_gene](auto copied_gene) {
                        return other_gene == copied_gene;
                    });
                });

                return *begin_lookup++;
            };

            std::generate(std::begin(sibling), std::begin(sibling) + low, lookup_unique);
            std::generate(std::begin(sibling) + high, std::end(sibling), lookup_unique);

            return sibling;
        };

        return { offspring(l, r), offspring(r, l) };
    }

    template <typename ChromoType>
    auto swap_mutation(ChromoType const& chromosome) -> ChromoType {
        std::size_t rnd1 = sto::random_between(0ul, std::size(chromosome)),
                    rnd2 = sto::random_between(0ul, std::size(chromosome));

        while (rnd2 == rnd1) {
            rnd2 = sto::random_between(0ul, std::size(chromosome));
        }

        auto mutated = chromosome;
        std::swap(mutated.at(rnd1), mutated.at(rnd2));

        return mutated;
    }

    template <typename ChromoType, typename GenOpFn>
    struct Genetic_Context {
        Population<ChromoType> m_current_population;

        std::vector<GenOpFn> m_phase1_ops;
        std::vector<GenOpFn> m_phase2_ops;
        std::vector<GenOpFn> m_phase3_ops;

        constexpr Genetic_Context(Population<ChromoType> initial_pop, std::vector<GenOpFn> phase1_ops,
                                    std::vector<GenOpFn> phase2_ops, std::vector<GenOpFn> phase3_ops)
            : m_current_population(initial_pop), m_phase1_ops(phase1_ops), m_phase2_ops(phase2_ops), m_phase3_ops(phase3_ops)
        {}

        auto evolve_once() -> Genetic_Context& {
            Population<ChromoType> phase1_population =
                std::accumulate(std::begin(m_phase1_ops), std::end(m_phase1_ops), m_current_population,
                            [](auto acc, auto op) {
                                    return op(acc);
                            });

            std::vector<Population<ChromoType>> phase2_pools;
            phase2_pools.reserve(std::size(m_phase2_ops));

            std::transform(std::begin(m_phase2_ops), std::end(m_phase2_ops), std::back_inserter(phase2_pools),
                            [&phase1_population](auto op) {
                               return op(phase1_population);
                            });

            Population<ChromoType> phase2_population =
                std::accumulate(std::begin(phase2_pools), std::end(phase2_pools), Population<ChromoType>{}, [](auto acc, auto cur) {
                    return acc + cur;
                });

            Population<ChromoType> phase3_population { phase2_population };

            while (std::size(phase3_population) < std::size(m_current_population)) {
                std::vector<Population<ChromoType>> phase3_pools;
                phase3_pools.reserve(std::size(m_phase3_ops));

                std::transform(std::begin(m_phase3_ops), std::end(m_phase3_ops), std::back_inserter(phase3_pools),
                               [&phase1_population](auto op) {
                                   return op(phase1_population);
                               });

                phase3_population =
                    std::accumulate(std::begin(phase3_pools), std::end(phase3_pools), phase3_population, [](auto acc, auto cur) {
                        return acc + cur;
                    });
            }

            m_current_population = phase3_population;

            return *this;
        }

        auto evolve_n(std::size_t n) -> Genetic_Context& {
            for (auto i = 0ul; i < n; ++i) { evolve_once(); }

            return *this;
        }

        template <typename Pred>
        auto evolve_until(Pred && pred) -> std::enable_if<std::is_same_v<std::invoke_result_t<Pred>, bool>, Genetic_Context&> {
            while (!pred(m_current_population)) { evolve_once(); };

            return *this;
        }
    };

    template <typename ChromoType, typename GenOpFn>
    Genetic_Context(Population<ChromoType> initial_pop, std::vector<GenOpFn> phase1_ops,
                        std::vector<GenOpFn> phase2_ops, std::vector<GenOpFn> phase3_ops)
                            -> Genetic_Context<ChromoType, GenOpFn>;

    template <typename ChromoType, typename GenOpFn>
    Genetic_Context(Population<ChromoType> initial_pop, std::initializer_list<GenOpFn> phase1_ops,
                    std::initializer_list<GenOpFn> phase2_ops, std::initializer_list<GenOpFn> phase3_ops)
    -> Genetic_Context<ChromoType, GenOpFn>;

    template <typename InputIt, typename ChromoType = std::remove_const_t<std::remove_reference_t<decltype(*std::declval<InputIt&>())>>>
    auto make_population(InputIt first, InputIt last) -> Population<ChromoType> {
        Population<ChromoType> pool;
        pool.reserve(std::distance(first, last));

        std::transform(first, last, std::back_inserter(pool),
                       [](auto chromosome) -> Organism<ChromoType> {
                           return { chromosome };
                       });

        return pool;
    }

    template <typename ChromoType>
    auto make_population(std::vector<ChromoType> container) -> Population<ChromoType> {
        return make_population(std::begin(container), std::end(container));
    }

    template <typename ChromoType>
    auto make_population(std::initializer_list<ChromoType> container) -> Population<ChromoType> {
        return make_population(std::begin(container), std::end(container));
    }
}

#endif // GENETICALGORITHM_HPP
