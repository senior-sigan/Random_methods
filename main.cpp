#include <QApplication>
#include "mainwindow.h"
#include "describing_statistic.h"
#include "statistical_criteria.h"
#include <iostream>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MainWindow w;
    w.show();
    return a.exec();
}
