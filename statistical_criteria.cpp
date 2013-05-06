#include "statistical_criteria.h"

statistical_criteria::statistical_criteria(QVector<qreal> in_rand_seq)
{
    this->rand_sequence = in_rand_seq;
}

void statistical_criteria::save(QString fname){
    QFile file(fname);
    file.open(QIODevice::WriteOnly | QIODevice::Text);
    QTextStream out(&file);
    out << "It works" << "\n";
    file.close();
}

qreal statistical_criteria::max(QVector<qreal> seq){
    qreal max = seq.at(0);
    for(int i=0;i<seq.size();i++) {
   //     std::cout<<"el="<< seq[i]<<std::endl;
        if (max < seq[i])
            max =  seq[i];
    }
    return max;
}

qreal statistical_criteria::min(QVector<qreal> seq){
    qreal min = seq.at(0);
    foreach (qreal el,seq) {
        if (el < min)
            min = el;
    }
    return min;
}

qreal statistical_criteria::phi(qreal z){

    qreal PI=3.141592653589;
    return (1/sqrt(2*PI))*exp(-pow(z,2)/2);

}

QVector<int> statistical_criteria::fill_the_array(int k,int h,int maxVal,int minVal){
    QVector<int> array;
    array.reserve(h);
    for(int i=0;i<h;i++){
        array.append(0);
    }
   foreach (qreal el,rand_sequence) {
        for (int i=0;i<maxVal-minVal;i+=k)
        {
            if ((el>=(minVal+i))&&(el<minVal+i+k)) {
                array[qFloor(i/k)]++;
                break;
            }
        }
    }
   return array;
}
QVector<qreal> statistical_criteria::find_the_sample_average_variance_sd(QVector<int>array,int k,int h,int minVal){

    QVector<qreal> answ;
    answ.reserve(2);
   QVector<qreal> relFreq;
   qreal sampleAverage=0;
   relFreq.reserve(h);
   for(int i=0;i<h;i++) {
       relFreq.append((qreal)array.at(i)/(qreal)this->rand_sequence.size());
       //  std::cout<<"relFreq.at(i)="<<relFreq.at(i)<<std::endl;
         sampleAverage+=((qreal)(minVal+k/2+i*k))*relFreq.at(i);
   }
   answ.append(sampleAverage);
   qreal sampleVariance=0;

   for(int i=0;i<h;i++) {
         sampleVariance+=pow((minVal+k/2+i*k-sampleAverage),2)*relFreq.at(i);
   }
  //  std::cout<<"sampleVariance="<<sampleVariance<<std::endl;
   qreal sampleSd=sqrt(sampleVariance);
   answ.append(sampleSd);
return answ;
}

QString statistical_criteria::chi_test(){
    QString res="";

    qreal maximum=max(rand_sequence);
    qreal minimum=min(rand_sequence);
    int minVal=qFloor(minimum);
    int maxVal=qCeil(maximum);
    int k=6;
    int h=qCeil((maxVal-minVal+1)/k);
    if ((qreal)(maxVal-minVal+1)/(qreal)k!=(qreal)h) h++;
    QVector<int> array;
    array.reserve(h);
    array=fill_the_array(k,h,maxVal,minVal);
/*
    std::cout<<"maxval="<<maxVal<<std::endl;
    std::cout<<"minval="<<minVal<<std::endl;
    std::cout<<"k="<<k<<std::endl;
    std::cout<<"h="<<h<<std::endl;
 */
  QVector<qreal> answ=find_the_sample_average_variance_sd(array,k,h,minVal);
  qreal sampleAverage=answ[0];
  qreal sampleSd=answ[1];
   QVector<qreal> zMas;
    zMas.reserve(h);
//std::cout<<"samplesd="<<sampleSd<<std::endl;
   for(int i=0;i<h;i++) {
         zMas.append((qreal)(minVal+k/2+i*k-sampleAverage)/sampleSd);
       // std::cout<<"zMas.at(i)="<<zMas.at(i)<<std::endl;
   }
   qreal sum=0;
   QVector<qreal> pMas;
    pMas.reserve(h);
    for(int i=0;i<h;i++) {
        pMas.append((qreal)k*(phi(zMas.at(i))/sampleSd));
        sum+=pMas[i];
        // std::cout<<"pMas.at(i)"<<pMas.at(i)<<std::endl;
    }
    std::cout<<sum<<std::endl;
   qreal khiRes=0;
   QVector<qreal> khi;
   khi.reserve(h);
   for(int i=0;i<h;i++) {
       khi.append((qreal)(pow((array[i]-this->rand_sequence.size()*pMas[i]),2)/(this->rand_sequence.size()*pMas[i])));
       khiRes+=khi[i];
   }
 std::cout<<"khiRes="<<khiRes<<std::endl;

 qreal table[31][13];
 QFile file("table_for_chi.txt");
 file.open(QIODevice::ReadOnly | QIODevice::Text);
 QTextStream file_to_read(&file);
 int i=0,j=0;
 while (!file_to_read.atEnd()){
     QString line=file_to_read.readLine();
     QStringList some_strings=line.split(" ");
     foreach(QString str,some_strings){
         table[i][j]=str.toDouble();
//         std::cout<<table[i][j]<<" ";
         j++;
     }
      //std::cout<<std::endl;
     i++;j=0;
  }
 file.close();

 int v=h-3;
 qreal alpha=0;
 if (v<=30){
     for (i=0;i<12;i++){
         if ((table[v][i]<=khiRes)&&(table[v][i+1]>=khiRes))
         {
             if (abs(table[v][i]-khiRes)<abs(table[v][i+1]-khiRes))
             alpha=table[0][i];
             else alpha=table[0][i+1];
             break;
         }
     }
 }
 else {
     qreal xMas[]={-2.33,-1.64,-0.674,0,0.674,1.64,2.33};
     qreal percent[]={0.99,0.95,0.75,0.5,0.25,0.05,0.01};
     qreal vych[7];
     for (i=0;i<7;i++){
         vych[i]=v+sqrt(2*v)*xMas[i]+(2/3)*pow(xMas[i],2)-2/3;
    std::cout<<" vych[i]="<< vych[i]<<std::endl;
     }
     for (i=0;i<6;i++){
         if ((vych[i]<=khiRes)&&(vych[i+1]>=khiRes))
         {
             if (abs(vych[i]-khiRes)<abs(vych[i+1]-khiRes))
             alpha=percent[i];
             else alpha=percent[i+1];
             break;
         }
     }

 }
  std::cout<<"alpha="<<alpha<<std::endl;
qreal alphaCr=0.03;
  if (alpha>alphaCr)std::cout<<"OKAAAAAY"<<std::endl;
  else std::cout<<"Everything is bad"<<std::endl;
    return res;
}


qreal statistical_criteria::distribution_curve(qreal ch,qreal table_raspr[41][10]){

    qreal ver;
    if (ch<5.0)
     ver=table_raspr[qFloor(ch*10)][(int)((ch*10-qFloor(ch*10))*10)];
    else ver=1;
    return ver;
}


QString statistical_criteria:: Kolmogorov_Smirnov_test(QVector<qreal> rand_seq){

    QVector<qreal> sequence_for_test=rand_seq;
    qreal maximum=max(sequence_for_test);
    qreal minimum=min(sequence_for_test);
    int minVal=qFloor(minimum);
    int maxVal=qCeil(maximum);
    int k=6;
    int h=qCeil((maxVal-minVal+1)/k);
    if ((qreal)(maxVal-minVal+1)/(qreal)k!=(qreal)h) h++;
  /*  std::cout<<"maxval="<<maxVal<<std::endl;
    std::cout<<"minval="<<minVal<<std::endl;
    std::cout<<"k="<<k<<std::endl;
    std::cout<<"h="<<h<<std::endl;*/
    QVector<int> array;
    array=fill_the_array(k,h,maxVal,minVal);
    QVector<qreal> answ=find_the_sample_average_variance_sd(array,k,h,minVal);
    qreal sampleAverage=answ[0];
    qreal sampleSd=answ[1];
    qSort(sequence_for_test);

    qreal table_raspr[41][10];
    QFile file("table_for_raspr.txt");
    file.open(QIODevice::ReadOnly | QIODevice::Text);
    QTextStream file_to_read(&file);
    int i=0,j=0;
    while (!file_to_read.atEnd()){
        QString line=file_to_read.readLine();
        QStringList some_strings=line.split(" ");
        foreach(QString str,some_strings){
            table_raspr[i][j]=str.toDouble();
        //    std::cout<<table_raspr[i][j]<<" ";
            j++;
        }
        // std::cout<<std::endl;
        i++;j=0;
     }
    file.close();
    QVector<qreal> plusMas,minusMas;
    int size=sequence_for_test.size();
    plusMas.reserve(size);
    minusMas.reserve(size);
    for (int i=0;i<size;i++){
        qreal value=(sequence_for_test[i]-sampleAverage)/sampleSd;
      //  std::cout<< "value="<<value<<std::endl;
        plusMas.append((qreal)(i+1)/(qreal)size-distribution_curve(value,table_raspr));
        minusMas.append(distribution_curve(value,table_raspr)-((qreal)i/(qreal)size));
    }
   // std::cout<< max(plusMas)<<std::endl;
  //  std::cout<< max(minusMas)<<std::endl;
    qreal plusK=sqrt(size)*max(plusMas);
    qreal minusK=sqrt(size)*max(minusMas);
    std::cout<<"plusK="<<plusK<<std::endl;
    std::cout<<"minusK="<<minusK<<std::endl;


    qreal table[13][7];
    QFile file2("table_for_kolm.txt");
    file2.open(QIODevice::ReadOnly | QIODevice::Text);
    QTextStream file_to_read2(&file2);
    i=0,j=0;
    while (!file_to_read2.atEnd()){
        QString line=file_to_read2.readLine();
        QStringList some_strings=line.split(" ");
        foreach(QString str,some_strings){
            table[i][j]=str.toDouble();
         //   std::cout<<table[i][j]<<" ";
            j++;
        }
     //    std::cout<<std::endl;
        i++;j=0;
     }
    file2.close();


    int n=size;
    qreal alpha1=0,alpha2=0;
    if (n<=12){
        for (i=0;i<6;i++){
            if ((table[n][i]<=plusK)&&(table[n][i+1]>=plusK))
            {
                if (abs(table[n][i]-plusK)<abs(table[n][i+1]-plusK))
                alpha1=table[0][i];
                else alpha1=table[0][i+1];
                break;
            }
            if ((table[n][i]<=minusK)&&(table[n][i+1]>=minusK))
            {
                if (abs(table[n][i]-minusK)<abs(table[n][i+1]-minusK))
                alpha2=table[0][i];
                else alpha2=table[0][i+1];
                break;
            }
        }
    }
    else {
        qreal yMas[]={0.07089,0.1601,0.3793,0.5887,0.8326,1.2239,1.5174};
        qreal percent[]={0.99,0.95,0.75,0.5,0.25,0.05,0.01};
        qreal vych[7];
        for (i=0;i<7;i++){
            vych[i]=yMas[i]-(qreal)(1.0/(6.0*sqrt(n)));
            std::cout<< vych[i]<<std::endl;
        }
        for (i=0;i<6;i++){
            if ((vych[i]<=plusK)&&(vych[i+1]>=plusK))
            {
                if (abs(vych[i]-plusK)<abs(vych[i+1]-plusK))
                alpha1=percent[i];
                else alpha1=percent[i+1];
                break;
            }

            if ((vych[i]<=minusK)&&(vych[i+1]>=minusK))
            {
                if (abs(vych[i]-plusK)<abs(vych[i+1]-plusK))
                alpha2=percent[i];
                else alpha2=percent[i+1];
                break;
            }
        }

    }

      qreal alphaCr=0.03;

      std::cout<<"alpha1="<<alpha1<<std::endl;

      if (alpha1>alphaCr)std::cout<<"OKAAAAAY"<<std::endl;
      else std::cout<<"Everything is bad"<<std::endl;
     std::cout<<"alpha2="<<alpha2<<std::endl;

     if (alpha2>alphaCr)std::cout<<"OKAAAAAY"<<std::endl;
     else std::cout<<"Everything is bad"<<std::endl;
   return "";
}

QString statistical_criteria::Subsequences_test(){

QVector<qreal> chetn,nech;

for(int i=0;i<this->rand_sequence.size();i++){
    if (i%2==0) chetn.append(this->rand_sequence.at(i));
    else nech.append(this->rand_sequence.at(i));
}
Kolmogorov_Smirnov_test(chetn);
Kolmogorov_Smirnov_test(nech);
return "";

}
