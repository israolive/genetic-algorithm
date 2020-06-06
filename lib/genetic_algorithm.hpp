#ifndef GENETICALGORITHM_HPP
#define GENETICALGORITHM_HPP

#include <algorithm>
#include <functional>
#include <vector>

namespace ai::gal {
    template <typename Chromosome>
    struct Organism {
        Chromosome c;
        double fitness {};
    };

    template <typename Chromosome>
    using Population = std::vector<Organism<Chromosome>>;

    auto rng() -> std::size_t {
        return static_cast<size_t>(std::rand());
    }

    template <typename Chromosome>
    struct Genetic_Operator {
        virtual auto operate(Population<Chromosome> const& pop) const -> Population<Chromosome> = 0;
    };

    template <typename Chromosome>
    auto operator|(Population<Chromosome> const& pop, Genetic_Operator<Chromosome> const& op) -> Population<Chromosome> {
        return op.operate(pop);
    }

    template <typename Chromosome>
    auto operator+(Population<Chromosome> & pool, Population<Chromosome> const& pop) -> Population<Chromosome>& {
        pool.reserve(pool.size() + pop.size());
        pool.insert(std::end(pool), std::begin(pop), std::end(pop));

        return pool;
    }

    namespace {
        template <typename Chromosome>
        inline auto random_binary_tournament(Population<Chromosome> const& population, double p, double best_chance = 0.8) -> std::size_t {
            std::size_t fst = rng() % population.size(),
                        snd = rng() % population.size();

            if (p < best_chance) {
                if (population.at(fst).fitness > population.at(snd).fitness) { return fst; } else { return snd; }
            } else {
                if (population.at(fst).fitness < population.at(snd).fitness) { return fst; } else { return snd; }
            }
        }
    }

    template <typename Chromosome>
    struct Tournament_Selector : public Genetic_Operator<Chromosome> {
        std::size_t m_count;
        double m_best_chance;

        Tournament_Selector(std::size_t count, double best_chance)
            : m_count(count), m_best_chance(best_chance) {}

        auto operate(Population<Chromosome> const& pop) const -> Population<Chromosome> override {
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

            Population<Chromosome> pool;
            pool.reserve(m_count);

            std::transform(std::begin(idxes), std::end(idxes), std::back_inserter(pool),
                           [&pop](auto idx) {
                               return pop.at(idx);
                           });

            return pool;
        }
    };

    template <typename Chromosome>
    struct Binary_Mating : public Genetic_Operator<Chromosome> {
        using Binary_Crossover_Fn = std::function<auto (Chromosome const&, Chromosome const&) -> std::pair<Chromosome, Chromosome>>;

        Binary_Crossover_Fn m_crossover;
        double m_probability;

        Binary_Mating(Binary_Crossover_Fn crossover_fn, double probability)
            : m_crossover(crossover_fn), m_probability(probability) {}

        auto operate(Population<Chromosome> const& pop) const -> Population<Chromosome> override {
            auto remainder = pop.size() % 2;

            Population<Chromosome> pool;
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

    template <typename Chromosome>
    struct Mutator : public Genetic_Operator<Chromosome> {
        using Mutator_Fn = std::function<auto (Chromosome const&) -> Chromosome>;

        Mutator_Fn m_mutator;
        double m_probability;

        Mutator(Mutator_Fn mutator_fn, double probability)
            : m_mutator(mutator_fn), m_probability(probability) {}

        auto operate(Population<Chromosome> const& pop) const -> Population<Chromosome> override {
            Population<Chromosome> pool;
            pool.reserve(pop.size());

            std::transform(std::begin(pop), std::end(pop), std::back_inserter(pool),
                           [mutator = m_mutator, p = m_probability](auto const& organism) {
                               if (rng() < p) {
                                   return Organism<Chromosome> { mutator(organism.c) };
                               } else {
                                   return organism;
                               }
                           });

            return pool;
        }
    };

}





#endif // GENETICALGORITHM_HPP
