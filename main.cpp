#include <QApplication>
#include "mainwindow.h"
#include "describing_statistic.h"
#include <iostream>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MainWindow w;
    w.show();
    QVector<qreal> lol;
    lol.append(1.0);
    lol.append(9.0);
    lol.append(5.0);
    lol.append(1.0);
    describing_statistic ds(lol);
    //std::cout<<ds.average()<<std::endl;
    ds.save("des_stat.txt");
    return a.exec();
}
