#include <vector>
#include "StitchPDF.h"
using namespace std;

extern "C" { 
    void estimatePDFstitch(double *sampleData, int *sampleLength, int *debug, 
                           int *outputLength, double *x, double *pdf, double *cdf,
                           int *resolution, int* smoothness, int* minLagrange){  
        OutputControl out;
        out.debug = debug[0];

        vector <double> samples;
        for (int i = 0; i < sampleLength[0]; i++) {
            samples.push_back(sampleData[i]);
        }

        InputParameters input;
        input.resolution = *resolution;
        input.minLagrange = *minLagrange;
        
        stitchPDF stitch = stitchPDF(samples, input);        
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
