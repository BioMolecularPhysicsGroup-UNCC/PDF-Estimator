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

class Block {
public:
    Block(vector <double> sample, int nPoints, int N, int blockNumber, bool debug, bool layerOne);
    virtual ~Block();       
    
    int blockNumber;
      
    vector <double> x;
    vector <double> pdf;
    vector <double> cdf;
        
    double xMin;
    double xMax; 
    
    double pdfPoint(double x);
    double cdfPoint(double x);    
    bool estimateBlock(double lowerBound, double upperBound);
    
    OutputControl out;   
   
private:
    
    vector <double> sample;    
    int nPoints;    
    int sampleLength;
    double normalize;
    bool layerOne;
    vector <double> xAll;
      
    WriteResults write;
    
    
};

#endif /* BLOCK_H */

