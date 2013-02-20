#include <QApplication>

#include <QHBoxLayout>

#include "Views/MainWindow.h"

#include "HCNGenerator.h"

#include "Views/WidgetHCN.h"

HCNGenerator *hcnGenerator;

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MainWindow w;
    w.resize(1024, 768);
    w.show();

    int mostCommonAlleleCountBins = 40;
    int distinctHaplotypeCountBins = 40;
    int sampleSize = 120;

//    MSModel model = MSModel::fromGrowthModel(25000, 50000, 0.1);
//    MSModel model = MSModel::fromGrowthModel(12500, 25000, 0.1);
//    MSModel model = MSModel::fromBottleneckModel(50000, 25000, 50000, 0.1, 0.2);
//    MSModel model = MSModel::fromGrowthModel(5000, 10000, 0.1);
    MSModel model = MSModel::fromBottleneckModel(10000, 2000, 10000, 0.1, 0.2);
    model.basepairs = 250000;

    hcnGenerator = new HCNGenerator(mostCommonAlleleCountBins, distinctHaplotypeCountBins, sampleSize, model);

    hcnGenerator->generate();

    w.widgetHCN->setResolution(mostCommonAlleleCountBins, distinctHaplotypeCountBins, &hcnGenerator->totalSamples);
    w.widgetHCN->setHCNData(hcnGenerator->hcn);

/*    for (int c = 0; c < mostCommonAlleleCountBins; ++c) {
        for (int d = 0; d < distinctHaplotypeCountBins; ++d) {
            qDebug("{%i %i %li}", c, d, hcnGenerator->hcn[c + d*mostCommonAlleleCountBins]);
        }
    }*/

    return a.exec();
}
