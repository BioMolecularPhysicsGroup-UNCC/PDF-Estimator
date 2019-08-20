/* 
 * File:   EstimatePDF.h
 * Author: jenny
 *
 * Created on February 2, 2019, 2:00 PM
 */

#ifndef ESTIMATEPDF_H
#define	ESTIMATEPDF_H

#include "mexAdapter.hpp"
#include "OutputControl.h"

using namespace matlab::data;

class MexFunction : public matlab::mex::Function {
    
private:
    std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr;        
    OutputControl out;
    
public:
    void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs);
};

#endif