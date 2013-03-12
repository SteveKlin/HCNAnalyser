#ifndef HCNGENERATOR_H
#define HCNGENERATOR_H

#include <algorithm>
#include <vector>

#include "Defines.h"

#include "MSSampleGenerator.h"

#define DBG(...)

// HCN summary results of a single simulated sample
class HCNSample {
private:
public:
    int distinctHaplotypeCount;
    int mostCommonAlleleFrequency;

    HCNSample(int iDistinctHaplotypeCount, int iMostCommonAlleleFrequency)
        : distinctHaplotypeCount(iDistinctHaplotypeCount), mostCommonAlleleFrequency(iMostCommonAlleleFrequency)
    {

    }
};

class HaplotypeComparator {
private:
    std::vector<int>& haplIndices;
    int usedSNPCount;
    int *snpSubsetIndices;
    char** list;
public:
    HaplotypeComparator(std::vector<int>& iHaplIndices, int* iSNPSubsetIndices, int iUsedSNPCount, char** iList) : haplIndices(iHaplIndices), snpSubsetIndices(iSNPSubsetIndices), usedSNPCount(iUsedSNPCount), list(iList) { }

    bool areEqual(int h, int g);

    bool operator() (int h, int g); // h lexicallylowerthan or equalto g
};

class HCNSampleGenerator {
private:
    int sampleSize;

    // Stats for snp count and accepted (satisfying the MAF constraint) snp counted
    long snpCountMin, snpCountMax;
    double snpCountTotal;
    int snpCountSampleCount;

    long accSNPCountMin, accSNPCountMax;
    double accSNPCountTotal;
    int accSNPCountSampleCount;
public:
    std::vector<HCNSample> samples;

    MSModel model;

    HCNSampleGenerator(int iSampleSize, MSModel iModel)
        : sampleSize(iSampleSize), model(iModel) { }

    ~HCNSampleGenerator() { }

    void generateSamples(int sampleCount);
};

class HCNGenerator {
private:
    int freqMostCommonBins, haplotypeCountBins;
    int haplotypeCountInitialStepSize;
    int sampleSize;
public:
    long totalSamples;
    long* hcn;

    MSModel model;

    HCNGenerator(int iFreqMostCommonBins, int iHaploTypeCountBins, int iSampleSize, MSModel iModel);

    void combineWithSamples(std::vector<HCNSample>& samples);

    void generate(int sampleCount = 10000);

    ~HCNGenerator() {
        delete [] hcn;
    }

    int getFreqMostCommonBin(int frequency) {
//        qDebug("gFMC(%i) %i %i", frequency, freqMostCommonBins, sampleSize);
        return (freqMostCommonBins*(frequency - 1))/sampleSize;
    }

    int getHaploTypeCountBin(int haplotypeCount) {
        return min(haplotypeCount/haplotypeCountInitialStepSize, haplotypeCountBins - 1);
    }

};

#endif // HCNGENERATOR_H
