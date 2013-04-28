#include "describing_statistic.h"

describing_statistic::describing_statistic(QVector<qreal> in_rand_seq){
    this->rand_sequence = in_rand_seq;
}
void describing_statistic::save(QString fname){
    QFile file(fname);
    file.open(QIODevice::WriteOnly | QIODevice::Text);
    QTextStream out(&file);
    out.setCodec("UTF8");
    QTextCodec *codec = QTextCodec::codecForName("UTF8");
    QTextCodec::setCodecForTr(codec);
    QTextCodec::setCodecForCStrings(codec);
    QTextCodec::setCodecForLocale(codec);

    out << "Count " << count()    << "\n"
        << "Sum " << sum() << "\n"
        << "Max " << max()      << "\n"
        << "Min " << min()      << "\n"
        << "Average " << average()  << "\n"
        << "Variance "<< variance() << "\n"
        << "Standard deviation " << standard_deviation() << "\n"
        << "Standard error " << standard_error() << "\n"
        << "Median " << median() << "\n"
        << "Moda " << moda_str() << "\n"
        << "Kurtosis " << kurtosis() << "\n"
        << "Skewness " << skewness() << "\n";
    file.close();
}
QString describing_statistic::moda_str(){
    QString res;
    foreach (qreal el, moda()){
        res.append(QString::number(el)+" ");
    }
    return res;
}

qreal describing_statistic::moment(int p){
    if (p<=0){
        throw "pow <= 0";
    }
    qreal aver = average();
    qreal mom = 0;
    int n = count();
    foreach (qreal y, this->rand_sequence){
        mom += pow(aver-y,p);
    }
    return mom/n;
}

qreal describing_statistic::variance(){
    return moment(2);
}
qreal describing_statistic::kurtosis(){
    return moment(4)/pow(variance(),2)-3;
}
QList<qreal> describing_statistic::moda(){
    QMap<qreal,int> modas;
    foreach (qreal el, this->rand_sequence){
        modas[el]++;
    }
    int maxim = modas.begin().value();
    foreach (qreal value, modas){
        if (value > maxim){
            maxim = value;
        }
    }
    QList<qreal> res;
    QMap<qreal,int>::iterator i;
    for (i=modas.begin();i != modas.end();++i){
        if (i.value() == maxim){
            res.append(i.key());
        }
    }
    return res;
}
qreal describing_statistic::median(qreal a){
    int n = count();
    qSort(rand_sequence);
    int k = a * (n -1);
    if ((k + 1) < a*n){
        return rand_sequence.at(k+1);
    }
    if ((k+1) == a*n){
        return (rand_sequence.at(k)+rand_sequence.at(k+1))/2;
    }
    return rand_sequence.at(k);
}
qreal describing_statistic::max(){
    qreal max = rand_sequence.at(0);
    foreach (qreal el,rand_sequence) {
        if (max < el)
            max = el;
    }
    return max;
}

qreal describing_statistic::min(){
    qreal min = rand_sequence.at(0);
    foreach (qreal el,rand_sequence) {
        if (min > el)
            min = el;
    }
    return min;
}
qreal describing_statistic::sum(){
    qreal s = 0;
    foreach (qreal el, rand_sequence) {
        s +=el;
    }
    return s;
}
qreal describing_statistic::average(){
    return sum()/ count();
}
qreal describing_statistic::skewness(){
    return moment(3)/pow(variance(),3.0/2.0);
}
qreal describing_statistic::standard_error(){
    return sqrt(variance()/count());
}

int describing_statistic::count(){
    return this->rand_sequence.size();
}
qreal describing_statistic::standard_deviation(){
    qreal var = variance();
    int n = count();
    qreal sd = 0;
    sd = n*var/(n-1);
    return sqrt(sd/n);
}
