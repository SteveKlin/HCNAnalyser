#include "WidgetHCN.h"

#include <QPainter>

#include "math.h"

WidgetHCN::WidgetHCN(QWidget *parent) : QWidget(parent) {
    hcnData = NULL;
    sampleCount = NULL;
    normalizedColors = true;
}

void WidgetHCN::setResolution(int newMostCommonAlleleCountBins, int newDistinctHaplotypeCountBins, long* newSampleCount) {
    if (hcnData) qDebug("Clear hcnData before changing hcn resolution!");
    distinctHaplotypeCountBins = newDistinctHaplotypeCountBins;
    mostCommonAlleleCountBins = newMostCommonAlleleCountBins;
    sampleCount = newSampleCount;
}

void WidgetHCN::setHCNData(long* newHCNData) {
    hcnData = newHCNData;
}

void WidgetHCN::paintEvent(QPaintEvent* paintEvent) {
    if (hcnData && (*sampleCount > 0)) {
        qDebug("w = %i, h = %i", width(), height());
        // Determine max value
        double normV = *sampleCount;
        if (normalizedColors) {
            long maxCount = 0;
            for (int d = 0; d < distinctHaplotypeCountBins; ++d) {
                for (int c = 0; c < mostCommonAlleleCountBins; ++c) {
                    long v = hcnData[c + d*mostCommonAlleleCountBins];
                    maxCount = (v > maxCount) ? v : maxCount;
                }
            }
            normV = (7/6.0)*maxCount;
        }


        QPainter painter(this);
        for (int d = 0; d < distinctHaplotypeCountBins; ++d) {
            for (int c = 0; c < mostCommonAlleleCountBins; ++c) {
                long v = hcnData[c + d*mostCommonAlleleCountBins];
                double vf = v/normV;
//                qDebug("%li %f", v, vf);
//                vf = vf/maxCountf;
                QColor kol = QColor::fromHslF(vf, 0.75, 0.66);
                painter.setBrush(QBrush(kol));
                painter.drawRect((c*width())/mostCommonAlleleCountBins, ((distinctHaplotypeCountBins - 1 - d)*height())/distinctHaplotypeCountBins,
                                 ((c + 1)*width())/mostCommonAlleleCountBins - (c*width())/mostCommonAlleleCountBins, ((distinctHaplotypeCountBins - 1 - d + 1)*height())/distinctHaplotypeCountBins - ((distinctHaplotypeCountBins - 1 - d)*height())/distinctHaplotypeCountBins
                                 );
            }
        }
    }
}
