#include "HCNGenerator.h"

// ### HaplotypeComparator #########################################

bool HaplotypeComparator::areEqual(int h, int g) {
//        int h = haplIndices[i];
//      int g = haplIndices[j];
    for (int c = 0; c < usedSNPCount; ++c) {
        if (list[h][snpSubsetIndices[c]] != list[g][snpSubsetIndices[c]]) return false;
    }
    return true;
}

bool HaplotypeComparator::operator() (int h, int g) { // h lexicallylowerthan or equalto g
//        int h = haplIndices[i];
//      int g = haplIndices[j];
    for (int c = 0; c < usedSNPCount; ++c) {
        if ((list[h][snpSubsetIndices[c]] == '1') && (list[g][snpSubsetIndices[c]] == '0')) return false;
        if ((list[h][snpSubsetIndices[c]] == '0') && (list[g][snpSubsetIndices[c]] == '1')) return true;
    }
    return true;
}

// ### HCNSampleGenerator ######################################

void HCNSampleGenerator::generateSamples(int sampleCount) {
    MSSampleGenerator sampleGenerator(sampleSize, model);

    int windowCount = 1;
    unsigned int snpSubsetSize = 20;
    unsigned int snpChoiceCount = 10;
    double mafLowerBound = 0.1;

    for (int sample = 0; sample < sampleCount; ++sample) {
//        qDebug("%i", sample);
        sampleGenerator.generateNextSample();

        double windowDefaultLength = 1.0/windowCount;
        int firstSNPOffset = 0;
        int lastSNPOffset = 0;
        for (int window = 0; window < windowCount; ++window) {
//              double windowStartPosition = window * windowDefaultLength;
//                double windowLength = windowDefaultLength;
            while ((lastSNPOffset < sampleGenerator.segsites) && (sampleGenerator.posit[lastSNPOffset] < (window + 1)*windowDefaultLength)) lastSNPOffset++;

            int snpCount = lastSNPOffset - firstSNPOffset + 1;

            // Make list of SNPs with MAF >= MAFLowerBound
            std::vector<int> acceptableSNPIndices(snpCount);
            int accIx = 0;
            for (int snp = 0; snp < snpCount; ++snp) {
                // Test if MAF >= MAFLowerBound
                int zeroC = 0, oneC = 0;
                for (int h = 0; h < sampleSize; ++h) {
                    if (sampleGenerator.list[h][firstSNPOffset + snp] == '0') { zeroC++; } else { oneC++; }
                }

                if (min(zeroC, oneC) >= mafLowerBound*sampleSize) {
                    acceptableSNPIndices[accIx++] = firstSNPOffset + snp;
                }
            }
            acceptableSNPIndices.resize(accIx);

//                qDebug("acceptable SNP count = %i", acceptableSNPIndices.size());

            unsigned int usedSNPCount = (acceptableSNPIndices.size() < snpSubsetSize) ? acceptableSNPIndices.size() : snpSubsetSize;
            for (unsigned int snpChoice = 0; snpChoice < snpChoiceCount; ++snpChoice) { // Different choices for reduction of Monte Carlo error of choice
                // Determine random choice of SNPs to use
                int snpSubsetIndices[usedSNPCount];

                if (usedSNPCount == acceptableSNPIndices.size()) {
                    for (unsigned int i = 0; i < usedSNPCount; ++i) snpSubsetIndices[i] = i;
                } else {
                    for (unsigned int i = 0; i < usedSNPCount; ++i) {
                        snpSubsetIndices[i] = (snpCount - i) * rand()/(RAND_MAX + 1.0); // Random choice between a numbered list of SNPs that have not been chosen already (in the original order)
                        // Adjust index to provide a direct index to an SNP that was not previously chosen
                        for (unsigned int j = 0; j < i; ++j) {
                            if (snpSubsetIndices[j] < snpSubsetIndices[i]) snpSubsetIndices[i]++;
                        }
                    }
                }

                // Adjust indices to map into original list
                for (unsigned int i = 0; i < usedSNPCount; ++i) snpSubsetIndices[i] = acceptableSNPIndices[snpSubsetIndices[i]];

                // Sort haplotypes so that we can efficiently determine window parameters
                std::vector<int> haplIndices(sampleSize);
                for (unsigned int s = 0; s < haplIndices.size(); ++s) haplIndices[s] = s;
                HaplotypeComparator haplotypeComparator(haplIndices, snpSubsetIndices, usedSNPCount, sampleGenerator.list);
                std::stable_sort(haplIndices.begin(), haplIndices.end(), haplotypeComparator);

                // Determine HCN coordinates for this window
                int maxHaplotypeCount = 1, currentHaplotypeCount = 1, distinctHaplotypeCount = 1;
                for (int h = 1; h < sampleSize; ++h) {
                    if (haplotypeComparator.areEqual(haplIndices[h - 1], haplIndices[h])) {
                        currentHaplotypeCount++;
                        if (currentHaplotypeCount > maxHaplotypeCount) maxHaplotypeCount = currentHaplotypeCount;
                    } else {
                        currentHaplotypeCount = 1;
                        distinctHaplotypeCount++;
                    }
                }

                samples.push_back(HCNSample(distinctHaplotypeCount, maxHaplotypeCount));

                if (usedSNPCount == acceptableSNPIndices.size()) break;
            }

            firstSNPOffset = lastSNPOffset + 1;
        }
    }

    sampleGenerator.printSegsiteStats();
}

// ### HCNGenerator ######################################################

HCNGenerator::HCNGenerator(int iFreqMostCommonBins, int iHaploTypeCountBins, int iSampleSize, MSModel iModel) {
    freqMostCommonBins = iFreqMostCommonBins;
    haplotypeCountBins = iHaploTypeCountBins;
    haplotypeCountInitialStepSize = 3;
    sampleSize = iSampleSize;
    model = iModel;

    hcn = new long[freqMostCommonBins * haplotypeCountBins];
    for (int f = 0; f < freqMostCommonBins; ++f) {
        for (int c = 0; c < haplotypeCountBins; ++c) {
            hcn[c*freqMostCommonBins + f] = 0;
        }
    }

    totalSamples = 0;
}

void HCNGenerator::combineWithSamples(std::vector<HCNSample>& samples) {
    for (unsigned int s = 0; s < samples.size(); ++s) {
        int c = getHaploTypeCountBin(samples[s].distinctHaplotypeCount);
        int f = getFreqMostCommonBin(samples[s].mostCommonAlleleFrequency);
   //     qDebug("[%i : %i :: %i : %i]", c, f, samples[s].distinctHaplotypeCount, samples[s].mostCommonAlleleFrequency);
        ensure(c >= 0);
        ensure(c < haplotypeCountBins);
        ensure(f >= 0);
        ensure(f < freqMostCommonBins);
        hcn[c*freqMostCommonBins + f]++;
        totalSamples++;
    }

/*        qDebug("(%i : %i)", freqMostCommonBins, haplotypeCountBins);

    for (int f = 0; f < freqMostCommonBins; ++f) {
        for (int c = 0; c < haplotypeCountBins; ++c) {
            qDebug("%i, %i: %i", f, c, hcn[c*freqMostCommonBins + f]);
        }
    }*/
}

void HCNGenerator::generate(int sampleCount) {
    HCNSampleGenerator sampleGenerator(sampleSize, model);

    sampleGenerator.generateSamples(sampleCount);

    combineWithSamples(sampleGenerator.samples);
}
