#include <vector>
#include "StitchPDF.h"
using namespace std;

extern "C" { 
    void estimatePDFstitch(double *sampleData, int *sampleLength, int *debug, 
                           int *outputLength, double *x, double *pdf, double *cdf){  
        OutputControl out;
        out.debug = debug[0];

        vector <double> samples;
        for (int i = 0; i < sampleLength[0]; i++) {
            samples.push_back(sampleData[i]);
        }        
        
        stitchPDF stitch = stitchPDF(samples);        
        stitch.out.debug = debug[0];
        stitch.process();
        vector <double> Vx = stitch.getX(outputLength[0]);
        vector <double> Vpdf = stitch.getPDF(outputLength[0]);  
        vector <double> Vcdf = stitch.getCDF(outputLength[0]);                     
        for (unsigned i = 0; i < Vx.size(); i++) {
            cdf[i] = Vcdf[i];
            pdf[i] = Vpdf[i];
            x[i] = Vx[i];
        }
        return;    
    }   
} 
