#include "view/MainWindow.hpp"

#include <QApplication>
#include <iostream>

#include "model/genetic_algorithm.hpp"

int main(int argc, char *argv[])
{
    using namespace ai::gal;

    QApplication a(argc, argv);
    Main_Window w;
    w.show();
    return a.exec();
}
