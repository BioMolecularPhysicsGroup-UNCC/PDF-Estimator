/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cppFiles/class.h to edit this template
 */

/* 
 * File:   Block.h
 * Author: jfarm
 *
 * Created on November 20, 2022, 11:45 AM
 */

#ifndef BLOCK_H
#define BLOCK_H

#include <vector>
#include <iostream>
#include "InputParameters.h"
#include "InputData.h"
#include "MinimizeScore.h"
#include "WriteResults.h"
#include "OutputControl.h"

using namespace std;

typedef struct LayerOptions {
    bool isLayerOne;
    bool isFarLeft;
    bool isFarRight;
} LayerOptions;

class Block {
public:
    Block(vector <double> sample, int nPoints, int N, int blockNumber, bool debug, LayerOptions opts);
    virtual ~Block();       
    
    int blockNumber;
      
    vector <double> x;
    vector <double> pdf;
    vector <double> cdf;
        
    double xMin;
    double xMax; 
    
    double pdfPoint(double x);
    double cdfPoint(double x);    
    double pdfPointHint(double x, int hint);
    double cdfPointHint(double x, int hint);
    bool estimateBlock(double lowerBound, double upperBound);
    
    OutputControl out;   
   
private:
    
    vector <double> sample;    
    int nPoints;    
    int sampleLength;
    double normalize;
    LayerOptions layerOpts;
    vector <double> xAll;
      
    WriteResults write;
    
    
};

inline double Block::pdfPointHint(double point, int i) {
    return normalize * (pdf[i - 1] + (point - xAll[i - 1]) * (pdf[i] - pdf[i - 1]) / (xAll[i] - xAll[i - 1]));
}

inline double Block::cdfPointHint(double point, int i) {
    return cdf[i - 1] + (point - xAll[i - 1]) * (cdf[i] - cdf[i - 1]) / (xAll[i] - xAll[i - 1]);
}

#endif /* BLOCK_H */

