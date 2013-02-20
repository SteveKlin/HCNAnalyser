#ifndef WIDGETHCN_H
#define WIDGETHCN_H

#include <QWidget>

class WidgetHCN : public QWidget {
    Q_OBJECT
private:
    int distinctHaplotypeCountBins, mostCommonAlleleCountBins;
    long *sampleCount;
    long* hcnData;
    bool normalizedColors;
public:
    explicit WidgetHCN(QWidget *parent = 0);
    
    void setResolution(int newMostCommonAlleleCountBins, int newDistinctHaplotypeCountBins, long* newSampleCount);
    void setHCNData(long* newHCNData);

    void paintEvent(QPaintEvent* paintEvent);
signals:
    
public slots:
    
};

#endif // WIDGETHCN_H
