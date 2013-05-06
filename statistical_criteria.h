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
    QVector<int> fill_the_array(int k,int h,int maxVal,int minVal);
    QVector<qreal> find_the_sample_average_variance_sd(QVector<int>array,int k,int h,int minVal);
    qreal phi(qreal z);
    qreal distribution_curve(qreal ch,qreal table_raspr[41][10]);
public:
    statistical_criteria(QVector<qreal> in_rand_seq);
    void save(QString fname);
    QString chi_test();
    QString Kolmogorov_Smirnov_test(QVector<qreal> rand_seq);
  //  QVector<qreal> sort_by_quadratic_selecion(QVector<qreal> mas, int size);
    QString Subsequences_test();
};

#endif // STATISTICAL_CRITERIA_H
