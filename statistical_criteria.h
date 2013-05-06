#ifndef STATISTICAL_CRITERIA_H
#define STATISTICAL_CRITERIA_H


#include <iostream>
#include <QtCore/qmath.h>
#include <QVector>
#include <QFile>
#include <QTextStream>
#include<QStringList>
#include <QtAlgorithms>

class statistical_criteria
{
    QVector<qreal> rand_sequence;
    qreal max(QVector<qreal> seq);
    qreal min(QVector<qreal> seq);
    QVector<long long> fill_the_array(long long k,long long h,long long maxVal,long long minVal);
    QVector<qreal> find_the_sample_average_variance_sd(QVector<long long>array,long long k,long long h,long long minVal);
    qreal phi(qreal z);
    qreal distribution_curve(qreal ch,qreal table_raspr[41][10]);
    QVector<QString> chi_test();
    QVector<QString> Kolmogorov_Smirnov_test(QVector<qreal> rand_seq);
    QVector<QString> Subsequences_test();

public:
    statistical_criteria(QVector<qreal> in_rand_seq);
    void save(QString fname,QString textToWrite);
    bool all_tests();
    void set_the_sequence(QVector<qreal> rand_seq);
    QVector<qreal> get_the_sequence();
};

#endif // STATISTICAL_CRITERIA_H
