#include "model/example/travelling_salesman.hpp"
#include "model/example/nqueens.hpp"

int main(int argc, char *argv[])
{
    std::srand(std::time(nullptr));

    std::puts("TSP"); travelling_salesman_problem();
    std::puts("NQueens"); nqueens_problem();
}