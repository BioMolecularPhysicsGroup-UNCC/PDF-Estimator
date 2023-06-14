
#include <vector>
#include "callPDF.h"

extern "C" { 
    void estimatePDF(double *sampleData, int *sampleLength, double *estimationPoints, double *estimationPDF, int *estimationLength, int *isSpecifyPoints,
                     double *low, double *high, int *isLow, int *isHigh, double * target,
                     int *lagrangeMin, int *lagrangeMax, int *debug, int *outlierCutoff, int *smooth,
                     int *points, int *fail, double *threshold, double *x, double *pdf, double *cdf, double *sqr, double *sqrSize, double *lagrange, double *r){    
        callPDF pd;
        pd.makeCall(sampleData, sampleLength[0], estimationPoints, estimationLength[0], isSpecifyPoints[0], low[0], high[0], isLow[0], isHigh[0], target[0], points[0], lagrangeMin[0], lagrangeMax[0], outlierCutoff[0], debug[0], smooth[0]);
        
        vector <double> Vx = pd.Vx;
        vector <double> Vpdf = pd.Vpdf;
        vector <double> VpdfPoints = pd.VpdfPoints;
        vector <double> Vcdf = pd.Vcdf;      
        vector <double> Vsqr = pd.Vsqr;  
        vector <double> Vlagrange = pd.Vlagrange;
        vector <double> Vr = pd.Vr;
        
        sqrSize[0] = pd.N;
        
        if (pd.solutionFailed) {
            fail[0] = 1;
        } else {            
            fail[0] = 0;
        }
        
        threshold[0] = pd.solutionThreshold;
        
        for (unsigned i = 0; i < Vx.size(); i++) {
            cdf[i] = Vcdf[i];
            pdf[i] = Vpdf[i];
            x[i] = Vx[i];
        }
        for (unsigned i = 0; i < VpdfPoints.size(); i++) {
            estimationPDF[i] = VpdfPoints[i];
        }
        for (unsigned i = 0; i < Vsqr.size(); i++) {
            sqr[i] = Vsqr[i];
        } 
        for (unsigned i = 0; i < Vr.size(); i++) {
            r[i] = Vr[i];
        }    

        for (unsigned i = 0; i < Vlagrange.size(); i++) {
            lagrange[i] = Vlagrange[i];
        }            
        return;    
    }   
} 
