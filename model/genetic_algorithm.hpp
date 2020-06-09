#ifndef GENETICALGORITHM_HPP
#define GENETICALGORITHM_HPP

#include <algorithm>
#include <numeric>
#include <functional>
#include <vector>

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

    auto rng() -> std::size_t {
        return static_cast<size_t>(std::rand());
    }

    template <typename ChromoType>
    using Genetic_Operator_Function = std::function<auto (Population<ChromoType> const&) -> Population<ChromoType>>;

    template <typename ChromoType>
    struct Genetic_Operator {

        virtual auto operate(Population<ChromoType> const& pop) const -> Population<ChromoType> = 0;

        virtual auto fun() -> Genetic_Operator_Function<ChromoType> {
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

    namespace {
        template <typename ChromoType>
        inline auto random_binary_tournament(Population<ChromoType> const& population, double p, double best_chance = 0.8) -> std::size_t {
            std::size_t fst = rng() % population.size(),
                        snd = rng() % population.size();

            if (p < best_chance) {
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
            std::vector<std::size_t> idxes (m_count);

            std::generate(std::begin(idxes), std::end(idxes),
                          [&idxes, &pop, chance = m_best_chance, last = std::begin(idxes)] () mutable {
                               std::size_t sel;

                               do {
                                   sel = random_binary_tournament(pop, rng(), chance);
                               } while (std::find(std::begin(idxes), last, sel) != last);

                               ++last;

                               return sel;
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

            for (std::size_t idx = 0; idx < pop.size() - remainder; idx += 2) {
                auto [first, second] = m_crossover(pop.at(idx).c, pop.at(idx + 1).c);

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
                               if (rng() < p) {
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

        Fitness_Applier(Fitness_Fn measurement)
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
    auto order_crossover(ChromoType const& l, ChromoType const& r) -> std::pair<ChromoType, ChromoType> {
        std::size_t rnd1 = rng() % std::size(l),
                    rnd2 = rng() % std::size(l);

        while (rnd2 == rnd1) {
            rnd2 = rng() % std::size(l);
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
        std::size_t rnd1 = rng() % std::size(chromosome),
                    rnd2 = rng() % std::size(chromosome);

        while (rnd2 == rnd1) {
            rnd2 = rng() % std::size(chromosome);
        }

        auto mutated = chromosome;
        std::swap(mutated.at(rnd1), mutated.at(rnd2));

        return mutated;
    }

    template <typename ChromoType, typename GenOpFn, std::size_t Phase1Count, std::size_t Phase2Count, std::size_t Phase3Count>
    struct Genetic_Context {
        Population<ChromoType> m_current_population;

        std::array<GenOpFn, Phase1Count> m_phase1_ops;
        std::array<GenOpFn, Phase2Count> m_phase2_ops;
        std::array<GenOpFn, Phase3Count> m_phase3_ops;

        Genetic_Context(Population<ChromoType> initial_pop,
                        std::array<GenOpFn, Phase1Count> phase1_ops,
                        std::array<GenOpFn, Phase2Count> phase2_ops,
                        std::array<GenOpFn, Phase3Count> phase3_ops)
            : m_current_population(initial_pop),
              m_phase1_ops(phase1_ops), m_phase2_ops(phase2_ops), m_phase3_ops(phase3_ops)
        {}

        auto evolve_once() -> Genetic_Context& {
            Population<ChromoType> phase1_population =
                std::accumulate(std::begin(m_phase1_ops), std::end(m_phase1_ops), m_current_population,
                            [](auto acc, auto op) {
                                    return op(acc);
                            });

            std::array<Population<ChromoType>, Phase2Count> phase2_pools;
            std::transform(std::begin(m_phase2_ops), std::end(m_phase2_ops), std::begin(phase2_pools),
                            [&phase1_population](auto op) {
                               return op(phase1_population);
                            });

            Population<ChromoType> phase2_population =
                std::accumulate(std::begin(phase2_pools), std::end(phase2_pools), Population<ChromoType>{}, [](auto acc, auto cur) {
                    return acc + cur;
                });

            Population<ChromoType> phase3_population { phase2_population };

            while (std::size(phase3_population) < std::size(m_current_population)) {
                std::array<Population<ChromoType>, Phase3Count> phase3_pools;
                std::transform(std::begin(m_phase3_ops), std::end(m_phase3_ops), std::begin(phase3_pools),
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
    };

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
