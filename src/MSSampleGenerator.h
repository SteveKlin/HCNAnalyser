#ifndef MSSAMPLEGENERATOR_H
#define MSSAMPLEGENERATOR_H

#include "../../ms/ms.h"
#include "ms-wrapper/mswrapper.h"

#include <stdio.h>
#include <QDebug>

#include <vector>
#include <string>

#include "math.h"


#include "Defines.h"

#include "stringTools.h"

extern int CI_MaxSites;

class DemographicEvent {
public:
    double time;
    double relativePopulationSize;

    DemographicEvent(double iTime, double iRelativePopulationSize) {
        time = iTime;
        relativePopulationSize = iRelativePopulationSize;
    }
};

class MSModel {
public:
    std::vector<DemographicEvent> demographicEvents;
    long currentPopulationSize;
    long basepairs;

    static MSModel fromGrowthModel(long Nmid, long Ncur, double tmid) {
        MSModel result;
        result.currentPopulationSize = Ncur;
        result.demographicEvents.push_back(DemographicEvent(tmid, Nmid/(double) Ncur));
        return result;
    }

    static MSModel fromBottleneckModel(long Nanc, long Nmid, long Ncur, double tmid, double tcur) {
        MSModel result;
        result.currentPopulationSize = Ncur;
        result.demographicEvents.push_back(DemographicEvent(tmid, Nmid/(double) Ncur));
        result.demographicEvents.push_back(DemographicEvent(tmid + tcur, Nanc/(double) Ncur));
        return result;
    }
};

class MSSampleGenerator {
private:
    double segfac;
    int count, ntbs, nseeds, maxSegSites;
    int numSamples;
    Parameters pars;

    int segsitesMin, segsitesMax;
    double segsitesTotal;
    double segsitesSampleCount;
//    struct params pars ;
public:
    int segsites;
    double probss, tmrca, ttot;
    double *posit;
    char **list;

    double thetaPopulationMutationRate; // = 4 * N_0 * basepairs * mu
    double muNeutralMutationRate;
    int basepairs;

    MSModel model;

    MSSampleGenerator(int iNumSamples, MSModel iModel, int segSitesIn = 0) {
        muNeutralMutationRate = 1E-9;
        model = iModel;
        thetaPopulationMutationRate = 4*model.currentPopulationSize*model.basepairs*muNeutralMutationRate;

        numSamples = iNumSamples;
        maxSegSites = (segSitesIn == 0) ? CI_MaxSites : segSitesIn;
        list = cmatrix(numSamples, maxSegSites + 1);
        posit = new double[maxSegSites];

        std::vector<std::string> arguments;
        arguments.push_back("ms");
        arguments.push_back(intToStr(numSamples));
        arguments.push_back(intToStr(1));
        arguments.push_back("-t");
        arguments.push_back(floatToStr(thetaPopulationMutationRate));

        for (unsigned int d = 0; d < model.demographicEvents.size(); ++d) {
            DemographicEvent& event = model.demographicEvents[d];
            arguments.push_back("-eN");
            arguments.push_back(floatToStr(event.time));
            arguments.push_back(floatToStr(event.relativePopulationSize));
        }

        char* args[arguments.size()];
        for (unsigned int a = 0; a < arguments.size(); ++a) {
            args[a] = new char[arguments[a].length() + 1];
            strcpy(args[a], arguments[a].c_str());
        }

        int howmany;
        getpars(arguments.size(), args, &howmany, &pars);
    }

    ~MSSampleGenerator() {
        delete [] posit;
        free_cmatrix(list, maxSegSites + 1);
    }

    void printSegsiteStats() {
        qDebug("segsitestats: min = %i, max = %i, avg = %f (%i)", segsitesMin, segsitesMax, segsitesTotal/(double) segsitesSampleCount);
    }

    void generateNextSample() {

/*                 if( pars.mp.theta > 0.0 ){
                    segfac = 1.0 ;
                    for(  i= pars.mp.segsitesin; i > 1; i--) segfac *= i ;
                 }  ?? */

        segsites = gensam(list, &probss, &tmrca, &ttot, pars, &posit);
//        qDebug("segsites = %i", segsites);
        if (segsitesSampleCount++ == 0) {
            segsitesMin = segsites;
            segsitesMax = segsites;
        } else {
            segsitesMin = min(segsites, segsitesMin);
            segsitesMax = max(segsites, segsitesMax);
        }
        segsitesTotal += segsites;
        /*
        if ( (segsites > 0 ) || ( pars.mp.theta > 0.0 ) ) {
           if( (pars.mp.segsitesin > 0 ) && ( pars.mp.theta > 0.0 )) fprintf(pf,"prob: %g\n", probss ) ;
           fprintf(pf,"segsites: %d\n",segsites);
           if( segsites > 0 )	fprintf(pf,"positions: ");
           for( i=0; i<segsites; i++) fprintf(pf,"%6.*lf ", pars.output_precision,posit[i] );
           fprintf(pf,"\n");
           if( segsites > 0 ) for(i=0;i<pars.cp.nsam; i++) { fprintf(pf,"%s\n", list[i] ); }
        }*/
    }
};

#endif // MSSAMPLEGENERATOR_H
