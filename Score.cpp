/* 
 * File:   Score.cpp
 * Author: jenny
 * 
 * Created on February 8, 2019, 8:33 PM
 */

#include "Score.h"


Score::Score(double confidenceTarget, double confidenceMin, double confidenceMax) {
    targetScore = getTargetScore(confidenceTarget);
    minimumScore = getTargetScore(confidenceMin);
    maximumScore = getTargetScore(confidenceMax);
}


Score::~Score() {
}


 
double Score::getTargetScore(double SURD) {
 
    vector<double>::iterator it;
    it = lower_bound (SURDs.begin(), SURDs.end(), SURD/100);
    unsigned index = it - SURDs.begin();
    
    if (index == SURDs.size()) {
        return scores[index - 1];
    } 
    if (index == 0) {
        return scores[0];
    } 
    
    double E1 = scores[index - 1];
    double E2 = scores[index];
    double P1 = SURDs[index - 1];
    double P2 = SURDs[index];
    double E = E1 + (SURD/100 - P1)*(E2 - E1)/(P2 - P1);
    return E;  
    
}

double Score::getConfidence(double score) {
 
    vector<double>::iterator it;
    it = lower_bound (scores.begin(), scores.end(), score);
    unsigned index = it - scores.begin();
    
    if (index == scores.size()) {
        return SURDs[index - 1];
    } 
    if (index == 0) {
        return SURDs[0];
    } 
    
    double E1 = scores[index - 1];
    double E2 = scores[index];
    double P1 = SURDs[index - 1];
    double P2 = SURDs[index];
    double P = P1 + (score - E1)*(P2 - P1)/(E2 - E1);
    return P * 100;  
    
}

