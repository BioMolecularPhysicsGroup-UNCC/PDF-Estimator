/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cppFiles/class.h to edit this template
 */

/* 
 * File:   stitchPDF.h
 * Author: jfarm
 *
 * Created on November 6, 2022, 10:25 AM
 */

#ifndef STITCHPDF_H
#define STITCHPDF_H

#include <vector>
#include <math.h>
#include <algorithm>
#include <iostream>
//#include <time.h>
#include "OutputControl.h"
#include "Block.h"

#define BUFFER (8)
#define RESOLUTION (50)

using namespace std;

class stitchPDF {
public:
    
    OutputControl out;
    
    stitchPDF(vector<double> sample, InputParameters input);
    stitchPDF(const stitchPDF& orig);
    virtual ~stitchPDF();    
        
    void process();
    
    vector <double> getX() {return x;};
    vector <double> getPDF() {return pdf;};
    vector <double> getCDF() {return cdf;};
    vector <double> getX(int nPoints);
    vector <double> getPDF(int nPoints);
    vector <double> getCDF(int nPoints);
    
    
    float getSortTime() {return duration;};
    
private:
    float duration;
    
    vector <double> x;
    vector <double> pdf;
    vector <double> cdf;
    
    const double thresholdCoeff = 0.5255;
    const double thresholdExp = 1.1;

    const double m0 =  1.09158290; 
    const double b0 = -0.17402986;
    const double m5 =  1.83306751; 
    const double b5 = -1.38564601;

    const int window = 10;
    const int maxLevel = 40;
    const int maxBlock = 100000;    
    
    int Ns;
    vector <double> sample;
    vector <int> partitions;
    
    vector <Block> blocks;

    InputParameters input;
    
    void branchNC();    
    
    void branch(int left, int right);
    void stagger();
    void estimate();
    void stitch();
    
    bool uniformSplit(int left, int right);
    inline double getRatio(const vector <double> &sample);
    inline double getThreshold(int n);
};

#endif /* STITCHPDF_H */

