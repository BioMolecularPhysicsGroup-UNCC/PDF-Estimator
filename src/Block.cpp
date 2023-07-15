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

Block::Block(vector <double> sample, int nPoints, int N, int blockNumber, bool debug, LayerOptions opts) {
    this->sample = sample;
    this->nPoints = nPoints;
    this->blockNumber = blockNumber;
    this->out.debug = debug;
    this->layerOpts = opts;
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
    input.SURDTarget = layerOpts.isLayerOne
        ? 50*(((double)sample.size() - 40.0)/99960.0) + 20 
        : 20;

    input.maxLagrange = (layerOpts.isLayerOne && !(layerOpts.isFarLeft || layerOpts.isFarRight))
        ? 8 
        : 250;

    input.smooth = 200;

    if (!layerOpts.isFarLeft) {
        input.lowerBound = lowerBound;
        input.lowerBoundSpecified = true;
    } else {
        input.lowerBoundSpecified = false;
    }

    if (!layerOpts.isFarRight) {
        input.upperBound = upperBound;
        input.upperBoundSpecified = true;
    } else {
        input.upperBoundSpecified = false;
    }

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
    long n = xAll.size();
    // 
    long i = min((long)(lower_bound(xAll.begin(), xAll.end(), point) - xAll.begin()), n - 1);  
    
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
    long n = xAll.size();
    double c; 
    if (point <= xAll[0]) {
        c = pdf[0];
    } else if (point >= xAll[n - 1]) {
        c = pdf[n - 1];
    } else {    
        long i = min((long)(lower_bound(xAll.begin(), xAll.end(), point) - xAll.begin()), n - 1);  

        c = pdf[i - 1] + (point - xAll[i - 1]) * (pdf[i] - pdf[i - 1]) / (xAll[i] - xAll[i - 1]);
    }

    return c * normalize;
}