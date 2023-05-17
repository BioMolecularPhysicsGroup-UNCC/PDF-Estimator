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

Block::Block(vector <double> sample, int nPoints, int N, int blockNumber) {
    this->sample = sample;
    this->nPoints = nPoints;
    this->blockNumber = blockNumber;
    sampleLength = sample.size();    
    normalize = 1.0 * sampleLength / N;
    estimateBlock();
}

Block::Block() {};

Block::~Block() {
}

void Block::estimateBlock() {              
    
    MinimizeScore minimumPDF = MinimizeScore();
    minimumPDF.out.debug = false;
    InputParameters input;
    input.SURDTarget = 20;
    input.maxLagrange = 100;
    input.smooth = 100;
    InputData data = InputData(input);
    
    data.setData(sample);     
    if (data.processData()) {      
        minimumPDF.minimize(input, data);
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
