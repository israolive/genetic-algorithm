#include <QApplication>
#include "view/MainWindow.hpp"

#include <iostream>
#include <array>
#include "model/genetic_algorithm.hpp"

int main(int argc, char *argv[])
{
    using namespace ai::gal;
    using namespace std::string_literals;

    //Chromosome<std::size_t> chr_a {1, 2, 3, 4, 5, 6, 7, 8, 9};
    //Chromosome<std::size_t> chr_b {5, 7, 4, 9, 1, 3, 6, 2, 8};
    std::string chr_a { "123456789" };
    std::string chr_b { "987654321" };

    auto pop = make_population({ chr_a, chr_b });

    Fitness_Applier<std::string> fitness([]([[maybe_unused]] std::string const& s) -> double {
        std::string const ideal = "541237698";

        double acc = 0;

        for (std::size_t i = 0; i < ideal.size(); i++) {
            if (s.at(i) == ideal.at(i)) {
                acc += 1;
            }
        }

        return acc;
    });

    Best_Selector<std::string> best_selector(1);

    Tournament_Selector<std::string> selector(2, 0.8);
    Binary_Mating<std::string> crossover(order_crossover<std::string>, 0.9);
    Mutator<std::string> mutator(swap_mutation<std::string>, 0.1);

    std::array phase1_ops { fitness.fun() };
    std::array phase2_ops { best_selector.fun() };
    std::array phase3_ops { static_cast<Genetic_Operator_Function<std::string>>([selector, crossover, mutator, fitness](auto const& pop) {
        return pop | selector | crossover | mutator | fitness;
    })};

    Genetic_Context<std::string, Genetic_Operator_Function<std::string>, 1, 1, 1> gc (pop, phase1_ops, phase2_ops, phase3_ops);

    auto evolved = gc.evolve_once().evolve_once();

    for (auto const& [ chromo, fitness ] : evolved.m_current_population) {
        std::cout << chromo << ' ' << fitness << std::endl;
    }

    School sc("aaasa");

    QApplication a(argc, argv);
    Main_Window w;
    w.show();
    return a.exec();
}
