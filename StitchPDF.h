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
#include <time.h>
#include "OutputControl.h"
#include "Block.h"


using namespace std;

class stitchPDF {
public:
    
    OutputControl out;
    
    stitchPDF(vector<double> sample);
    stitchPDF(const stitchPDF& orig);
    virtual ~stitchPDF();    
        
    void process();
    
    vector <double> getX() {return x;};
    vector <double> getPDF() {return pdf;};
    vector <double> getCDF() {return cdf;};
    
    float getSortTime() {return duration;};
    
private:
    float duration;
    
    vector <double> x;
    vector <double> pdf;
    vector <double> cdf;
    
    const double thresholdCoeff = 0.5255;
    const double thresholdExp = 1.1;
    const int window = 10;
    const int maxLevel = 40;
    const int maxBlock = 100000;    
    
    int Ns;
    vector <double> sample;
    vector <int> partitions;
    
    vector <Block> blocks;
    
    void branchNC();    
    
    void branch(int left, int right);
    void stagger();
    void estimate();
    void stitch();
    
    int uniformSplit(vector <double> sample);
    double getRatio(vector <double> sample, int start, int stop);
};

#endif /* STITCHPDF_H */
