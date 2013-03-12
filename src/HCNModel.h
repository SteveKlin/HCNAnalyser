#ifndef HCNMODEL_H
#define HCNMODEL_H

#include <vector>

#include "HCNGenerator.h"
#include "MSSampleGenerator.h"

class HCNModel {
public:
    HCNModel();

    virtual void generateHCNs() { }
};

class HCNGrowthModelParameters {
public:
    double tMid;
};

template <class ModelParameters>
class ModelData {
public:
    ModelParameters parameters;
    long* hcn;

    ModelData() : hcn(NULL) { }
    ~ModelData() { /*delete hcn;*/ }

    void storeHCN(long* newHCN, int mostCommonAlleleCountBins, int distinctHaplotypeCountBins) {
  //      hcn = newHCN;
        hcn = new long[mostCommonAlleleCountBins*distinctHaplotypeCountBins];
        memcpy(hcn, newHCN, sizeof(hcn[0]) * mostCommonAlleleCountBins * distinctHaplotypeCountBins);
        adjust(mostCommonAlleleCountBins, distinctHaplotypeCountBins);
    }

    void adjust(int mostCommonAlleleCountBins, int distinctHaplotypeCountBins) {
        int binCount = mostCommonAlleleCountBins * distinctHaplotypeCountBins;
        for (int q = 0; q < binCount; ++q) hcn[q]++;
    }

    double getLikelihoodFor(long* otherHCN, int mostCommonAlleleCountBins, int distinctHaplotypeCountBins) {
        int binCount = mostCommonAlleleCountBins * distinctHaplotypeCountBins;
        long sampleCount = 0;
        for (int q = 0; q < binCount; ++q) {
            sampleCount += hcn[q];
        }

        double probFactor = 1.0/sampleCount;
        double likelihood = 0.0;
        for (int q = 0; q < binCount; ++q) {
//            qDebug("[%li]", hcn[q]);
            likelihood += otherHCN[q] * log(probFactor * hcn[q]);
        }
        return likelihood;
    }
};

class HCNGrowthModel : public HCNModel {
private:
    std::vector<ModelData<HCNGrowthModelParameters> > data;
    int mostCommonAlleleCountBins, distinctHaplotypeCountBins, sampleSize;
    int NMid, NCur, basepairs;
    int resolutionTCur;
    double startTMid, endTMid;
public:
    HCNGrowthModel() {
        NMid = 1000;
        NCur = 3000;
        basepairs = 250000;
        mostCommonAlleleCountBins = 16;
        distinctHaplotypeCountBins = 16;
        sampleSize = 120;
    }

    double getTMidForModel(int index) {
        return startTMid + index*(endTMid - startTMid)/resolutionTCur;
    }

    void generateHCNs(double iStartTMid, double iEndTMid, int iResolutionTCur) {
        startTMid = iStartTMid;
        endTMid = iEndTMid;
        resolutionTCur = iResolutionTCur;
        for (int q = 0; q < resolutionTCur; ++q) {
            double tMid = getTMidForModel(q);
            qDebug("Generating model (%i/%i) tMid = %f", q + 1, resolutionTCur, tMid);
            MSModel msModel = MSModel::fromGrowthModel(NMid, NCur, tMid);
            msModel.basepairs = basepairs;

            HCNGenerator hcng(mostCommonAlleleCountBins, distinctHaplotypeCountBins, sampleSize, msModel);
            hcng.generate(2500);

            ModelData<HCNGrowthModelParameters> d;
            d.parameters.tMid = tMid;
            d.storeHCN(hcng.hcn, mostCommonAlleleCountBins, distinctHaplotypeCountBins);
            data.push_back(d);
//            d.hcn = NULL;
        }
    }

    void estimateParameters(double tMid) {
        MSModel msModel = MSModel::fromGrowthModel(NMid, NCur, tMid);
        msModel.basepairs = basepairs;

        HCNGenerator hcng(mostCommonAlleleCountBins, distinctHaplotypeCountBins, sampleSize, msModel);
        hcng.generate(5000);

        double optLikelihood;
        int optLIndex;
        for (int i = 0; i < data.size(); ++i) {
            double hcnLikelihood = data[i].getLikelihoodFor(hcng.hcn, mostCommonAlleleCountBins, distinctHaplotypeCountBins);
            qDebug("Testing model (%i/%i): i = %i LL = %f", i + 1, data.size(), i, hcnLikelihood);
            if ((i == 0) || (hcnLikelihood > optLikelihood)) {
                optLIndex = i;
                optLikelihood = hcnLikelihood;
            }
        }
        qDebug("Opt = %i LL = %f  E[tMid] = %f", optLIndex, optLikelihood, getTMidForModel(optLIndex));
    }
};

#endif // HCNMODEL_H
