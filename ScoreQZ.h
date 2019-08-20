/* 
 * File:   ScoreQZ.hpp
 * Author: jenny
 *
 * Created on February 8, 2019, 8:33 PM
 */

#ifndef SCOREQZ_HPP
#define	SCOREQZ_HPP

#include "Score.h"

using namespace std;


class ScoreQZ : public Score{
public:
    ScoreQZ(double confidenceTarget, double confidenceMin, double confidenceMax); 
    virtual ~ScoreQZ();
    virtual vector <int> getIndices (int N, int p);
    virtual vector <int> setIndices (int N, int p);
    double calculateScore(double r[], int N, int p);
    void getValues();   
    
};

#endif	/* SCORE_HPP */

