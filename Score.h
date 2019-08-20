/* 
 * File:   Score.hpp
 * Author: jenny
 *
 * Created on February 8, 2019, 8:33 PM
 */

#ifndef SCORE_HPP
#define	SCORE_HPP

#include <fstream>
#include <vector>
#include <iostream>
#include <math.h>


#include <algorithm>

using namespace std;

class Score {
public:
    double targetScore;
    double minimumScore;
    double maximumScore;
    
    Score() {};
    Score(const Score& orig) {};
    Score(double confidenceTarget, double confidenceMin, double confidenceMax); 
    virtual ~Score();
    
    virtual double calculateScore(double r[], int N, int p) {return 0;}
    virtual vector <int>  setIndices (int N, int p) {}
    virtual vector <int> getIndices (int N, int p) {return indices;}
    void setSpacing(int N, int P) {}
    
    double getLikelihood() {return likelihood;};
    void setVarianceMin(bool qzVar) {minimizeVariance = qzVar;}
    double getConfidence(double score);    
    double SURD;    
    double QZVariance = 0;
    
protected:
    vector <double> scores;
    vector <double> SURDs;
    double likelihood;
    vector <int> indices;
    
    int spacing;
    
    bool minimizeVariance = true;
    
    double getTargetScore(double SURD);  
    const double PI = 3.141592;    
    void getValues();
    
};

#endif	/* SCORE_HPP */

