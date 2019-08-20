/* 
 * File:   PDFMainMatlab.cpp
 * Author: jenny
 * 
 * Created on December 3, 2018, 8:32 AM
 */

#include "EstimatePDF.h"
#include "callPDFMatlab.h" 

    
void MexFunction::operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {
        
    
    int nParameters = inputs.size();    
    matlabPtr = getEngine();
    out.debug = true;
    
    matlab::data::ArrayFactory factory;
    if (nParameters < 1) {        
        out.displayError(matlabPtr, "Data sample required");
    }
    if (nParameters > 2) {        
        out.displayError(matlabPtr, "Only two input arguements allowed:  data sample and parameter structure");
    }
            
    TypedArray<double> doubleArray = std::move(inputs[0]);
    int sampleSize = doubleArray.getNumberOfElements();
    double sampleData[sampleSize];
    int i = 0;
    for (auto& elem : doubleArray) {
        sampleData[i++] = elem;
    }                  
     
    
    callPDF *pd = new callPDF(matlabPtr);
    
    if (nParameters > 1) {
        StructArray const matlabStructArray = inputs[1];
        auto fields = matlabStructArray.getFieldNames();
        std::vector<std::string> fieldNames(fields.begin(), fields.end());
            
        int count = 0;
        Array structField;
        for (std::vector<string>::iterator iter = fieldNames.begin(); iter != fieldNames.end(); ++iter) {
            string field = *iter;          
            if (strcmp(field.c_str(), "SURDtarget") == 0) {                
                structField = matlabStructArray[0][fieldNames[count++]];                
                pd->setSURDtarget(structField[0]);
            } else if (strcmp(field.c_str(), "SURDmin") == 0) {                
                structField = matlabStructArray[0][fieldNames[count++]];                
                pd->setSURDmin(structField[0]);
            } else if (strcmp(field.c_str(), "SURDmax") == 0) {                
                structField = matlabStructArray[0][fieldNames[count++]];                
                pd->setSURDmax(structField[0]);        
            } else if (strcmp(field.c_str(), "scoreType") == 0) {                
                structField = matlabStructArray[0][fieldNames[count++]];                
                pd->setScoreType((CharArray(structField)).toAscii());
            } else if (strcmp(field.c_str(), "LagrangeMin") == 0) {                
                structField = matlabStructArray[0][fieldNames[count++]];                
                pd->setLagrangeMin(structField[0]);
            } else if (strcmp(field.c_str(), "LagrangeMax") == 0) {                
                structField = matlabStructArray[0][fieldNames[count++]];                
                pd->setLagrangeMax(structField[0]);
            } else if (strcmp(field.c_str(), "integrationPoints") == 0) {                
                structField = matlabStructArray[0][fieldNames[count++]];                
                pd->setPoints(structField[0]);
            } else if (strcmp(field.c_str(), "lowBound") == 0) {                
                structField = matlabStructArray[0][fieldNames[count++]];                
                pd->setLow(structField[0]);
            } else if (strcmp(field.c_str(), "highBound") == 0) {                
                structField = matlabStructArray[0][fieldNames[count++]];                
                pd->setHigh(structField[0]);
            } else if (strcmp(field.c_str(), "outlierCutoff") == 0) {
                structField = matlabStructArray[0][fieldNames[count++]];
                pd->setOutlierCutoff(structField[0]);
            } else if (strcmp(field.c_str(), "fuzz") == 0) {                
                structField = matlabStructArray[0][fieldNames[count++]];                
                pd->setFuzz(structField[0]);
            } else if (strcmp(field.c_str(), "partition") == 0) {                
                structField = matlabStructArray[0][fieldNames[count++]];                
                pd->setPartitionSize(structField[0]);
            } else if (strcmp(field.c_str(), "debug") == 0) {                
                structField = matlabStructArray[0][fieldNames[count++]];                
                pd->setDebug(structField[0]);  
            } else if (strcmp(field.c_str(), "minVariance") == 0) {                
                structField = matlabStructArray[0][fieldNames[count++]];                
                pd->setVariance(structField[0]);
            } else if (strcmp(field.c_str(), "adaptiveDx") == 0) {
                structField = matlabStructArray[0][fieldNames[count++]];
                pd->setAdaptiveDx(structField[0]);
            } else {
                string unknown = "Unknown parameter: " + field;
                out.displayError(matlabPtr, unknown);
            }                
        }  
    }                
    
       
    pd->makeCall(sampleSize, sampleData);            
    
    vector <int> failed;
    failed.push_back(pd->solutionFailed);
    ArrayDimensions sz = {failed.size(), 1, 1};
    TypedArray<int> MsolutionFailed = factory.createArray(sz, failed.begin(), failed.end());
    outputs[0] = MsolutionFailed;  
           
    if (!pd->solutionFailed) {                
        sz = {pd->Vx.size(), 1, 1};        
        TypedArray<double> Mx = factory.createArray(sz, pd->Vx.begin(), pd->Vx.end());
        outputs[1] = Mx;
        
        sz = {pd->Vpdf.size(), 1, 1};
        TypedArray<double> Mpdf = factory.createArray(sz, pd->Vpdf.begin(), pd->Vpdf.end());
        outputs[2] = Mpdf;
        
        sz = {pd->Vcdf.size(), 1, 1};
        TypedArray<double> Mcdf = factory.createArray(sz, pd->Vcdf.begin(), pd->Vcdf.end());
        outputs[3] = Mcdf;
        
        sz = {pd->Vsqr.size(), 1, 1};
        TypedArray<double> Msqr = factory.createArray(sz, pd->Vsqr.begin(), pd->Vsqr.end());
        outputs[4] = Msqr;
        
        sz = {pd->Vlagrange.size(), 1, 1};
        TypedArray<double> Mlagrange = factory.createArray(sz, pd->Vlagrange.begin(), pd->Vlagrange.end());
        outputs[5] = Mlagrange;                     
        
        vector <double> threshold;
        threshold.push_back(pd->thresholdScore);
        sz = {threshold.size(), 1, 1};
        TypedArray<double> Mthreshold = factory.createArray(sz, threshold.begin(), threshold.end());
        outputs[6] = Mthreshold;  
    
        vector <double> confidence;
        confidence.push_back(pd->confidenceScore);
        sz = {confidence.size(), 1, 1};
        TypedArray<double> Mconfidence = factory.createArray(sz, confidence.begin(), confidence.end());
        outputs[7] = Mconfidence;  
        
        sz = {pd->Vr.size(), 1, 1};
        TypedArray<double> Mr = factory.createArray(sz, pd->Vr.begin(), pd->Vr.end());
        outputs[8] = Mr;
    
    } 
    delete pd;             
}
    
    
    
