/* 
 * File:   OutputControl.h
 * Author: jenny
 *
 * Created on October 31, 2018, 6:57 PM
 */


#ifndef OUTPUTCONTROL_H
#define	OUTPUTCONTROL_H

//#define outputCommandLine
#define outputMatlab
//#define outputR


#include <string>

#ifdef outputMatlab
#include "cppmex/mexMatlabEngine.hpp"
#include "MatlabDataArray/ArrayFactory.hpp"
#endif

#ifdef outputR
#include "R_ext/Print.h"
#endif

#ifdef outputCommandLine
#include <iostream>
#endif

using namespace std;
class OutputControl {
public:
    OutputControl();
    OutputControl(const OutputControl& orig);
    virtual ~OutputControl();    
    bool debug;
    
#ifdef outputMatlab    
    std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr; 
    
    void displayError(std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr, std::string errorMessage);
    void displayOnMATLAB(std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr, std::string message);  
    void setPtr(std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr); 
#endif
    
    
void print(string output);
void print(string output, int value);
void print(string output, double value);
void error(string output);
void error(string output, int value);
void error(string output, double value);       

};

#endif	/* OUTPUTCONTROL_H */


