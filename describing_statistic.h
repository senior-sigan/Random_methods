#ifndef DESCRIBING_STATISTIC_H
#define DESCRIBING_STATISTIC_H

#include <QCoreApplication>
#include <QVector>
#include <math.h>
#include <QMap>
#include <QFile>
#include <QTextStream>
#include <QTextCodec>
#include <QTranslator>

class describing_statistic
{
    QVector<qreal> rand_sequence;

public:
    describing_statistic(QVector<qreal> in_rand_seq);//Конструктор на входе получает вектор случайных чисел, а не файл!!
    qreal sum();//сумма последовательности
    int count();//длина последовательности
    qreal max();
    qreal min();
    qreal moment(int p);//p-ый момент
    qreal kurtosis();//эксцесс
    QList<qreal> moda();//мода
    QString moda_str();
    qreal median(qreal a=0.5);//медиана
    qreal variance();//дисперсия
    qreal skewness();//ассиметричность
    qreal average();//среднее
    qreal standard_deviation();//стандартное отклонение
    qreal standard_error();
    void save(QString fname); //Сохранить всю статистику в файл
};

#endif // DESCRIBING_STATISTIC_H
