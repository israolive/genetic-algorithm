#include "MainWindow.hpp"
#include "ui_MainWindow.h"

Main_Window::Main_Window(QWidget *parent)
    : QMainWindow(parent)
      , ui(new Ui::Main_Window)
{
    ui->setupUi(this);
}

Main_Window::~Main_Window()
{
    delete ui;
}

