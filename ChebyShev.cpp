/* 
 * PDF Estimator:  A non-parametric probability density estimation tool based on maximum entropy
 * File:   ChebyShev.cpp
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

#include "ChebyShev.h"

ChebyShev::ChebyShev() {
   
}

ChebyShev::ChebyShev(const ChebyShev& orig) {
}

ChebyShev::~ChebyShev() {
}

void ChebyShev::initialize(double dzLocal[], int sizeLocal) {
    this->size = sizeLocal;
    this->dz = dzLocal;
    vector <double> zeroT;
    vector <double> oneT;
    vector <double> twoT;

    for (int z = 0; z < size; z++) {
        double x = -1 + dz[z]*2;
        zeroT.push_back(1);
        oneT.push_back(x);
        twoT.push_back(2*x*oneT[z] - 1);
    }
    termsT.push_back(zeroT);
    termsT.push_back(oneT);
    termsT.push_back(twoT);
}

double* ChebyShev::getTerms(unsigned mode) {   
    if (termsT.size() <= mode) {
        vector <double> test = addMode(mode);
        return &test[0];
    } else {
        return &termsT[mode][0];
    }
}

vector < vector < double > > ChebyShev::getAllTerms(unsigned mode) {
    for (unsigned i = 0; i < mode; i++) {
        if (termsT.size() <= i) {
            addMode(i);
        }
    }
    return termsT;
}


vector <double> ChebyShev::addMode(int mode) {        
         
        vector <double> T = termsT.at(mode-1);
        vector <double> Tprev = termsT.at(mode-2);
        vector <double> Tnext; 
        double x = 0;       
               
        for (int z = 0; z < size; z++) {
            x = -1 + 2*dz[z];
            Tnext.push_back(2*x*T[z] - Tprev[z]);
        }
        termsT.push_back(Tnext);
        return Tnext;
    }
    