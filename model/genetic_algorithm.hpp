#ifndef GENETICALGORITHM_HPP
#define GENETICALGORITHM_HPP

#include <algorithm>
#include <functional>
#include <vector>

namespace ai::gal {
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
    struct Genetic_Operator {
        virtual auto operate(Population<ChromoType> const& pop) const -> Population<ChromoType> = 0;
    };

    template <typename ChromoType>
    auto operator|(Population<ChromoType> const& pop, Genetic_Operator<ChromoType> const& op) -> Population<ChromoType> {
        return op.operate(pop);
    }

    template <typename ChromoType>
    auto operator+(Population<ChromoType> & pool, Population<ChromoType> const& pop) -> Population<ChromoType>& {
        pool.reserve(pool.size() + pop.size());
        pool.insert(std::end(pool), std::begin(pop), std::end(pop));

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

            std::transform(std::begin(idxes), std::end(idxes), std::begin(idxes),
                           [&idxes, &pop, chance = m_best_chance, last = std::begin(idxes)] ([[maybe_unused]] auto holder) mutable {
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

    template <typename Allele>
    using Chromosome = std::vector<Allele>;

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
}





#endif // GENETICALGORITHM_HPP
