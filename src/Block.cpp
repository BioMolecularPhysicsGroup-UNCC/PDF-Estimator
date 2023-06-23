/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cppFiles/class.cc to edit this template
 */

/* 
 * File:   Block.cpp
 * Author: jfarm
 * 
 * Created on November 20, 2022, 11:45 AM
 */

#include "Block.h"

Block::Block(vector <double> sample, int nPoints, int N, int blockNumber, bool debug, bool layerOne) {
    this->sample = sample;
    this->nPoints = nPoints;
    this->blockNumber = blockNumber;
    this->out.debug = debug;
    this->layerOne = layerOne;
    out.print("in block");
    sampleLength = sample.size();    
    normalize = 1.0 * sampleLength / N;
    //estimateBlock();
}

Block::~Block() {
}

bool Block::estimateBlock(double lowerBound, double upperBound) {              
    
    MinimizeScore minimumPDF = MinimizeScore();
    minimumPDF.out.debug = false;//out.debug;
    InputParameters input;

    // dynamic SURD Target
    input.SURDTarget = layerOne ? 50*(((double)sample.size() - 40.0)/99960.0) + 20 : 20;

    input.maxLagrange = layerOne ? 8 : 250;
    input.smooth = 200;

    input.lowerBound = lowerBound;
    input.upperBound = upperBound;

    input.upperBoundSpecified = true;
    input.lowerBoundSpecified = true;
    InputData data = InputData(input);
    
    out.print("sample size ", (int) sample.size());
    
    data.setData(sample);     
    if (data.processData()) {      
        bool minimized_failed = minimumPDF.minimize(input, data);
        write.createSolution(input, data, minimumPDF);        
        xAll = write.x;
        pdf = write.PDF;  
        cdf = write.CDF;  
        nPoints = xAll.size();
        
        vector <double> gridPoints;
        gridPoints.reserve(nPoints);
        double dx = 1.0 / (nPoints - 1);
        double sum = 0;//dx;
        for (int i = 0; i < nPoints - 2; i++) {
            gridPoints.push_back(sum);
            sum += dx;
        }      
        gridPoints.push_back(1.0);      
      
        x = write.interpolateGrid(cdf, write.x, gridPoints);  
        xMin = write.min;
        xMax = write.max;

        return !minimized_failed;
       
/*        ostringstream blockName; 
        blockName << blockNumber;
        string filename;       
        
        filename = "blocktest_" + blockName.str() + ".txt";
        vector <double> pdfInt;
        for (int i = 0; i < x.size(); i++) {
            pdfInt.push_back(pdfPoint(x[i]));
        }
 
        write.writeColumns(filename, x, pdfInt, x.size());
 */
  
        
    } else {
        return false;
    }
}
 
double Block::cdfPoint(double point) {
    int n = xAll.size();
    int i = 0;   
    
    while (xAll[i] < point) {
        if (++i == (n - 1)){
            break;
        } 
    }
    double c;
    if ((i == 0) || (i >= (n - 1))) {
        if (i == 0) {
            c = 0;
        } else {
            c = 1;
        }
    } else {
        c = cdf[i - 1] + (point - xAll[i - 1]) * (cdf[i] - cdf[i - 1]) / (xAll[i] - xAll[i - 1]);
    }
    return c;
}

double Block::pdfPoint(double point) {
    int n = xAll.size();
    double c;
    int i = 0;   
    if (point <= xAll[i]) {
        c = pdf[i];
    } else {
        if (point >= xAll[n - 1]) {
            c = pdf[n - 1];
        } else {    
            while (xAll[i] < point) {
                if (++i == (n - 1)){
                    break;
                } 
            }
            c = pdf[i - 1] + (point - xAll[i - 1]) * (pdf[i] - pdf[i - 1]) / (xAll[i] - xAll[i - 1]);
        }
    }
    return c * normalize;
}
