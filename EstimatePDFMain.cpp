/* 
 * PDF Estimator:  A non-parametric probability density estimation tool based on maximum entropy
 * File:   Estimate.cpp
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

#include <time.h>

#include <cstdlib>
#include <vector>
#include <stdio.h>
#include <fstream>
#include "InputParameters.h"
#include "Variable.h"
#include "JointProbability.h"
#include "InputData.h"
#include "ScoreQZ.h"
#include "MinimizeScore.h"
#include "WriteResults.h"
#include "StitchPDF.h"
using namespace std;

int main(int argc, char** argv) {               
    bool stitch = true;
    double readTime;                                                  
    readTime = time(NULL);
    OutputControl out;   
    out.debug = true;
    InputParameters input;   
    if (!input.userInput(argc, argv)) {
        return false;
    }   
   
    ifstream fin;
    string line;
    fin.open((input.inputPath + input.inputFile).c_str());
    if(!fin.is_open()){
        out.error("Failed to open data file " + input.inputFile);
        return false;
    }
    vector <double> sample;
    string temp;
    int nVariables;   
    vector < vector < double > > rawData;     
    bool firstTime = true;
    int nColumns = 0;
    double test;
    while (getline(fin, line)) {
        nColumns = 0;
        stringstream checkLine(line);
        while (getline(checkLine, temp, ' ')) {
            nColumns++;
            test = atof(temp.c_str());
            if (!firstTime) {
                if (nColumns > nVariables) {
                    int row = rawData.size();
                    out.error("too many input columns in line  ", row);
                    return false;
                }
                sample[nColumns - 1] = test;
            } else {                
                sample.push_back(test);
            }
            
        }
        if (firstTime) {
            firstTime = false;
            nVariables = nColumns;
        }
        if (nColumns > nVariables) {
            int row = rawData.size();
            out.error("too few input columns in line  ", row);
            return false;
        } 
        rawData.push_back(sample);
    }
    if (rawData.size() < 10) {
        out.error("Sample must contain at least 10 data points");
        return false;        
    }   
    fin.close();       
            
    vector <Variable> variables;
    if (!stitch) {
        variables.reserve(nVariables);
    }
        
    int row = rawData.size();    
    for (int v = 0; v < nVariables; v++) {
        vector <double> samples;
        for (int i = 0; i < row; i++) {
            samples.push_back(rawData[i][v]);
        }               
        readTime = time(NULL) - readTime;        
        if (!stitch) {
            ostringstream vString; 
            vString << v; 
            Variable variable = Variable(input, samples, vString.str(), false);
            variables.push_back(variable);
        } else {        
            stitchPDF stitch = stitchPDF(samples);                          
            double stitchTime = time(NULL);
            stitch.process();  
            stitchTime = time(NULL) - stitchTime;   
            WriteResults write;
            write.writeColumns((input.inputPath + "out_" + input.inputFile).c_str(), stitch.getX(), stitch.getPDF(), stitch.getX().size());
        }
    }              

    if (!stitch) {
        int nPoints = (int) (200 + rawData.size()/200.0);
        if (nPoints > 1500) nPoints = 1500;
        double allTime;                                                  
        allTime = time(NULL);
        JointProbability jp = JointProbability(variables, row, nPoints);
   
        jp.calculate();
        vector <double> x = jp.getRange();
        vector <double> pdf = jp.getJP();
        allTime = time(NULL) - allTime;        
        out.print("total time", allTime);    
//    WriteResults write;
//    write.writeValuesAppend((input.inputPath + "out_" + input.inputFile).c_str(), stitchTime, allTime);
//    write.writeColumns((input.inputPath + "out_" + input.inputFile).c_str(), x, pdf, x.size());
    }
    return 0;
}
