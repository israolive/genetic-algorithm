#include "view/MainWindow.hpp"

#include <QApplication>
#include <iostream>

#include "model/genetic_algorithm.hpp"

int main(int argc, char *argv[])
{
    using namespace ai::gal;

    Chromosome<std::size_t> chr_a {1, 2, 3, 4, 5, 6, 7, 8, 9};
    Chromosome<std::size_t> chr_b {5, 7, 4, 9, 1, 3, 6, 2, 8};

    [[maybe_unused]] auto [ sib, _ ] = order_crossover(chr_a, chr_b);

    for (auto gene : sib) {
        std::cout << gene << ' ';
    }

    std::cout << std::endl;

    QApplication a(argc, argv);
    Main_Window w;
    w.show();
    return a.exec();
}
