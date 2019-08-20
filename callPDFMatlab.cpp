/* 
 * File:   callPDFMatlab.cpp
 * Author: jenny
 * 
 * Created on December 3, 2018, 8:36 AM
 */


#include "callPDFMatlab.h"


callPDF::callPDF(std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr) {
    input = new InputParameters;
    this->matlabPtr = matlabPtr;
    out.setPtr(matlabPtr);
    out.debug = debug;    
}

callPDF::callPDF(const callPDF& orig) {
}

callPDF::~callPDF() {    
    delete input;
}

void callPDF::setLow(double low) {
    input->lowerBoundSpecified = true;
    input->lowerBound = low;
}

void callPDF::setHigh(double high) {
    input->upperBoundSpecified = true;
    input->upperBound = high;
}

void callPDF::setSURDtarget(int surd) {
    input->SURDTarget = surd;
}

void callPDF::setSURDmin(int surd) {
    input->SURDMinimum = surd;
}

void callPDF::setSURDmax(int surd) {
    input->SURDMaximum = surd;
}

void callPDF::setScoreType(string type) {
    input->scoreType = type;
}

void callPDF::setPoints(int points) {
    input->integrationPoints = points - 1;
}

void callPDF::setFuzz(bool fuzz) {
    input->fuzz = fuzz;
}

void callPDF::setDebug(bool debug) {
    this->debug = debug;
}

void callPDF::setPartitionSize(int partition) {
    input->initPartitionSize = partition;
}

void callPDF::setLagrangeMin(int lagrange) {
    input->minLagrange = lagrange;
}

void callPDF::setLagrangeMax(int lagrange) {
    input->maxLagrange = lagrange;
}

void callPDF::setOutlierCutoff(double cutoff) {
    input->outlierCutoff = cutoff;
}

void callPDF::setAdaptiveDx(bool adaptive) {
    input->adaptive = adaptive;
}

void callPDF::makeCall(int sampleLength, double* sampleData) {
   
    
    out.debug = true;
    Score* score;
    if (input->scoreType == "AD") {
        ScoreAD *scoreAD = new ScoreAD(input->SURDTarget, input->SURDMinimum, input->SURDMaximum);
        score = scoreAD;
    } else if (input->scoreType == "LL") {               
        ScoreLL *scoreLL = new ScoreLL(input->SURDTarget, input->SURDMinimum, input->SURDMaximum);
        score = scoreLL;
    } else {
        ScoreQZ *scoreQZ = new ScoreQZ(input->SURDTarget, input->SURDMinimum, input->SURDMaximum);
        score = scoreQZ;
    }  
    score->setVarianceMin(minVariance);
    MinimizeScore *minimumPDF = new MinimizeScore();
    InputData *data = new InputData(*input);  
    
    data->out.setPtr(matlabPtr);
    data->out.debug = debug;
    minimumPDF->out.setPtr(matlabPtr);
    minimumPDF->out.debug = debug;
    input->out.setPtr(matlabPtr);
    input->out.debug = debug;
    vector <double> inputData;
    for (int i = 0; i < sampleLength; i++) {
        inputData.push_back(sampleData[i]);
    }
    input->writeHeader = false;
    input->writeFile = false;
    if (input->initPartitionSize == 0) {
        input->initPartitionSize = sampleLength;
    }
    data->setData(inputData);     
    if (data->processData()) {
        solutionFailed = minimumPDF->minimize(input, *data, *score);
        if (!solutionFailed) {            
            WriteResults write;        
            write.out.debug = debug;
            write.out.setPtr(matlabPtr);
            write.createSolution(input, data, minimumPDF, score);  
//            Vlagrange = minimumPDF->getLagrange();        
            write.createQQ(minimumPDF->trialRandom, minimumPDF->N);  
            thresholdScore = minimumPDF->bestThreshold;
            confidenceScore = score->getConfidence(minimumPDF->bestThreshold);
            double * r = minimumPDF->trialRandom;
            for (int i = 0; i < minimumPDF->N; i++) {
                Vr.push_back(r[i]); 
            }
            Vlagrange = write.L;
            Vsqr = write.SQR;
            Vcdf = write.CDF;
            Vpdf = write.PDF;
            Vx   = write.x;     
        }
    } else {
        solutionFailed = true;
    }
    
    delete data;
    delete score;
    delete minimumPDF;
    
}

