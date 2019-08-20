/* 
 * PDF Estimator:  A non-parametric probability density estimation tool based on maximum entropy
 * File:   inputData.cpp
 * Copyright (C) 2018
 * Jenny Farmer jfarmer6@uncc.edu
 * Donald Jacobs djacobs1@uncc.edu
 * 
 * This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published 
 * by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in 
 * the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
 * PURPOSE.  See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with 
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <numeric>

#include "InputData.h"
#include "WriteResults.h"


InputData::InputData(const InputParameters& input) {
    this->input = input;
        
    N = 0;
    nPoints = 0;
    minimumRaw = 0;
    maximumRaw = 0;
    minimumCalc = 0;
    maximumCalc = 0;
    nPointsAdjust = 0;
    
    nRightOutliers = 0;
    nLeftOutliers = 0;    
    leftOutliers = false;
    rightOutliers = false;
    
}

InputData::InputData(const InputData& orig) {
}

InputData::~InputData() {    
    delete [] doubleInverse;
    delete [] transformedZeroOne;
    delete [] inverse;
    delete [] dz;
    delete [] xUntransform;    
}

bool InputData::readData() {
    
    ifstream fin;
    string line;
    fin.open((input.inputPath + input.inputFile).c_str());

    if(!fin.is_open()){
        out.error("Failed to open data file " + input.inputFile);
        return false;
    }
	
    int temp = 0;
    while (getline(fin, line)) {
        double test = atof(line.c_str());
        if (test == 0) {
            test = 0;
        }
        temp++;
        rawData.push_back(test);
    }
    if (rawData.size() == 0) {
        out.error("No data in " + input.inputFile);
        return false;        
    }
   
    fin.close();   
    return processData();
}
    
void InputData::setData(vector<double> data) {
    rawData.resize(data.size());
    rawData = data;
}

bool InputData::processData() {   
    nPoints = input.integrationPoints;
    if (nPoints == -1) {
        nPoints = (int) (200 + rawData.size()/200.0);
        if (nPoints > 1500) nPoints = 1500;
    }
    
    sort(rawData.begin(), rawData.end());
    minimumRaw = rawData[0];
    maximumRaw = rawData[rawData.size() - 1];
    if (minimumRaw == maximumRaw) {        
        out.error("All input data has the same value ", minimumRaw);
        return false;        
    }
    
    int nValues = rawData.size();
    if (input.upperBoundSpecified) {  
            maximumCalc = input.upperBound;
    } else {
        double max = rawData[nValues - 1];
        maximumCalc = max + (max - rawData[rawData.size() - 5]);
    }
    if (input.lowerBoundSpecified) {  
        minimumCalc = input.lowerBound;
    } else {
        double min = rawData[0];
        minimumCalc = min + (min - rawData[4]);
    }
    
    if (input.outlierCutoff > 0) {        
        identifyOutliers();
    }
    
    if (!transformData()) {
        return false;
    }
    setAdaptiveDz();  
    cheby.initialize(doubleInverse, 2*nPointsAdjust-1);
    return true;
}


  void InputData::identifyOutliers() {
         
        double q1 = 0;
        double q3 = 0;      
        
        int nValues = rawData.size();    
        
        int middle = (int) (nValues/2);
        int quarter = (int) (middle/2);
        
        if (nValues%2 == 0) {
            if (middle%2 == 0) {
                q1 = (rawData[(int)quarter - 1] + rawData[(int)quarter])/2;
                q3 = (rawData[(int)quarter + (int)middle - 1] + rawData[(int)quarter + (int)middle])/2;
            } else {
                q1 = rawData[(int)quarter];
                q3 = rawData[(int)quarter + (int)middle];
            }
        } else {
            if (middle%2 == 0) {
                q1 = (rawData[(int)quarter - 1] + rawData[(int)quarter])/2;
                q3 = (rawData[(int)quarter + (int)middle] + rawData[(int)quarter + (int)middle + 1])/2;
            } else {
                q1 = rawData[(int)quarter];
                q3 = rawData[(int)quarter + (int) middle] + 1;
            }
        }               
        double iqr = input.outlierCutoff*(q3 - q1);
        double leftOutlier = q1 - iqr;
        double rightOutlier = q3 + iqr;                   
        
        if (maximumCalc > rightOutlier) {
            maximumCalc = rightOutlier;
            rightOutliers = true;                
        }
        if (minimumCalc < leftOutlier) {
            minimumCalc = leftOutlier;
            leftOutliers = true;                
        }
  }
    


bool InputData::transformData() {                    
    
    
    int nValues = rawData.size();     
            
    for (vector<double>::iterator iter = rawData.begin(); iter != rawData.end(); ++iter) {
        if (*iter >= minimumCalc) {
            if (*iter <= maximumCalc) {
                tempData.push_back(*iter);
            } else {
                nRightOutliers++;
            }
        } else {
            nLeftOutliers++;
        }            
    }
        
    nValues = tempData.size(); 
    if (nValues == 0) {
        out.error("No data within specified boundaries");
        return false;
    }
    
    
    transformedData.clear();      
    transformedData.reserve(nValues);
    transformedZeroOne = new double[nValues];
    int count = 0;
    for (vector<double>::iterator iter = tempData.begin(); iter != tempData.end(); ++iter) {
        transformedData.push_back((2*(*iter) - maximumCalc - minimumCalc)/(maximumCalc - minimumCalc));
        transformedZeroOne[count] = (transformedData[count] + 1)/2.0;
        count++;    
    }
    return true;
}


void InputData::setAdaptiveDz() {
   
    vector <double> dzVector;   
    N = transformedData.size();
     
    double dzMax = 2.0/(nPoints - 1);
        
    int skip = (int) (N/(nPoints - 1));
    if (skip==0) skip = 1;
                        
    double last = -1.0;
    double next;
        
    for (int b = skip; b <= (N + skip); b+=skip) {
        if (b >= (N)) {
            next = transformedData[N-1];
        }
        else {
            next = transformedData[b];
        }
        double test = next - last;
        double difference = fabs(test);
        if (difference > dzMax) {
            double steps =  difference/dzMax;
            int iSteps = (int) steps;
            for (int k = 0; k < (iSteps + 1); k++) {
                dzVector.push_back(difference/(iSteps + 1));
            }
        }
        else {             
            dzVector.push_back(difference); 
        }            
        last = next;
    }            
        
    int dzSize = dzVector.size();
    inverse = new double[dzSize + 1];
    inverse[0] = 0;
    for (int j = 1; j <= dzSize; j++) {
        inverse[j] = inverse[j-1] + dzVector[j-1]/2.0;
    }

    double difference = 1.0 - inverse[dzSize];
    if (difference > dzMax) {
        double steps =  (difference)/dzMax;
        int iSteps = (int) steps;
        for (int k = 0; k < 2*(iSteps + 1); k++) {
            dzVector.push_back(difference/(iSteps + 1));
        }        
    }
    else {
        dzVector.push_back(difference); 
    }
    dzSize = dzVector.size();
    dz = new double[dzSize];
    delete [] inverse;    
    inverse = new double[dzSize];
    doubleInverse = new double[2*dzSize - 1];
    inverse[0] = dzVector[0]/2.0;
    dz[0] = dzVector[0]/2.0;
    int count = 0;
    for (int j = 1; j < dzSize; j++) {
        dz[j] = dzVector[j]/2.0;
        inverse[j] = inverse[j-1] + dzVector[j]/2.0;
        doubleInverse[count] = inverse[j-1];
        doubleInverse[count+1] = (inverse[j-1] + inverse[j])/2.0;
        count += 2;
    }        
    doubleInverse[count] = (inverse[dzSize-1] + 1.0)/2.0;
    
    xUntransform = new double[2*dzSize - 1];
        
    for (int i=0; i < 2*dzSize - 1; i++) {
        xUntransform[i] = doubleInverse[i]*2.0 - 1;
    }
    for (int i=0; i < 2*dzSize - 1; i++) {
        xUntransform[i] = (maximumCalc - minimumCalc)*xUntransform[i] + minimumCalc + maximumCalc;
        xUntransform[i] /= 2;
    }
        
    
    nPointsAdjust = dzSize;        
}
    
    
