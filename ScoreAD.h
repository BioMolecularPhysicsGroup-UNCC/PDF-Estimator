/* 
 * File:   ScoreAD.hpp
 * Author: jenny
 *
 * Created on February 4, 2019, 3:31 PM
 */

#ifndef SCOREAD_HPP
#define	SCOREAD_HPP
#include "Score.h"

using namespace std;

class ScoreAD : public Score {
public:
    ScoreAD(double confidenceTarget, double confidenceMin, double confidenceMax); 
    virtual ~ScoreAD();
    virtual vector <int> getIndices (int N, int p);
    virtual vector <int> setIndices (int N, int p);
    double calculateScore(double r[], int N, int p);
    void getValues();   
};

#endif	/* SCOREAD_HPP */

