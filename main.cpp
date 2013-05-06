#include <QApplication>
#include "mainwindow.h"
#include "describing_statistic.h"
#include "statistical_criteria.h"
#include "gauss_generator.h"
#include <iostream>


QVector<qreal> generate_the_sequence(){
    generate_gauss_seq("test.txt", 10, 10000, 69069, 5,((unsigned int)-1));
    QVector<qreal> sequence;
    QFile file2("test.txt");
    file2.open(QIODevice::ReadOnly | QIODevice::Text);
    QTextStream file_to_read2(&file2);
    int i=0;
    while (!file_to_read2.atEnd()){
        QString line;
        file_to_read2>>line;
        QStringList some_strings=line.split(" ");
        foreach(QString str,some_strings){
           sequence.append(str.toDouble());
           i++;
        }
     }
    file2.close();
    sequence.remove(i-1);
    return sequence;
}

void explore_the_sequence(QVector<qreal> sequence){
    statistical_criteria st(sequence);
    bool result=st.all_tests();
    while (result==false){
        QVector<qreal> sequence=generate_the_sequence();
        st.set_the_sequence(sequence);
        result=st.all_tests();
    }
    describing_statistic ds(st.get_the_sequence());
    ds.save("ds.txt");
}
int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MainWindow w;
    w.show();
    QVector<qreal> sequence=generate_the_sequence();
    explore_the_sequence(sequence);
    return a.exec();
}
