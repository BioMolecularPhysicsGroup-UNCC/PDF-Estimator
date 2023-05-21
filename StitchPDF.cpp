#include <omp.h>
#include <iostream>

/* 
 * File:   stitchPDF.cpp
 * Author: jfarm
 * 
 * Created on November 6, 2022, 10:25 AM
 */

#include "StitchPDF.h"

stitchPDF::stitchPDF(vector<double> sample) {
    
    Ns = sample.size();    
    sort(sample.begin(), sample.end());       
    this->sample = sample;
}

void stitchPDF::process() {
    clock_t algorithmTime;                                                  
    algorithmTime = clock();
    branch(0, Ns - 1);  
    algorithmTime = clock() - algorithmTime;  
    duration = ((float) algorithmTime)/CLOCKS_PER_SEC;
    algorithmTime = clock();
    stagger(); 
    estimate();   
    stitch(); 
    algorithmTime = clock() - algorithmTime;  
    duration = ((float) algorithmTime)/CLOCKS_PER_SEC;
}

stitchPDF::stitchPDF(const stitchPDF& orig) {
}

stitchPDF::~stitchPDF() {
}


void stitchPDF::branch(int left, int right) {
    vector <double> range(sample.begin() + left, sample.begin() + right + 1);
    int partition = uniformSplit(range);
    if (partition > 0) {        
        if ((right - (partition + left)) > window) {
            branch(left, partition + left);
        }
        partition += left;
        partitions.push_back(partition);
        if ((right - partition) > window) {
            branch(partition, right);
        }
    }
}


int stitchPDF::uniformSplit(vector <double> sample) {    
    int n = sample.size();    
    double threshold = thresholdCoeff * pow(n, thresholdExp);    
    int partition = (int) ceil(n / 2);
    
    vector <double> dx;
    dx.reserve(n - 1);
    for (int i = 0; i < (n - 1); i++) {
        dx.push_back(sample[i + 1] - sample[i]);
    }  
    
    double ratioLeft = getRatio(dx, 0, partition - 1);
    double ratioRight = getRatio(dx, partition - 1, n - 1);     
    double average = (ratioLeft + ratioRight) / 2;
    if (average < threshold) {
        partition = -1;             
    }    
    return partition;
}

double stitchPDF::getRatio(vector <double> dxAll, int start, int stop) {
    int n = stop - start + 1;    
    vector <double> dxSort;
    dxSort.reserve(dxAll.size());
    dxSort = dxAll;
    sort(dxSort.begin() + start, dxSort.begin() + stop);
    double dxMin = 0;
    double dxMax = 0;
    for (int i = 0; i < window; i++) {
        dxMin += dxSort[i + start];
        dxMax += dxSort[stop - i - 1];
    }       
    return dxMax/dxMin;    
}


void stitchPDF::stagger() {
    partitions.push_back(Ns - 1);
    partitions.insert(partitions.begin(), 0);
    int nPartitions = partitions.size();
    for (int p = 0; p < nPartitions - 1; p++) {
        partitions.push_back((int) ceil((partitions[p] + partitions[p+1]) / 2));
    }
    sort(partitions.begin(), partitions.end());
    
    nPartitions = partitions.size();
    
    int steps;
    for (int p = 0; p < nPartitions - 1; p++) {  
        double blockSize = partitions[p + 1] - partitions[p];
        if (blockSize > maxBlock) {
            steps = (int) ceil(blockSize / maxBlock);
            double add = blockSize / (steps);
            for (int iSteps = 0; iSteps < (steps - 1); iSteps++) {                
                partitions.push_back(partitions[p] + add * (iSteps + 1));
            }
        }
    }
    sort(partitions.begin(), partitions.end());
}

void stitchPDF::estimate() {
    int nPartitions = partitions.size();    
  
    blocks.reserve(nPartitions - 2);
    blocks.resize(nPartitions - 2);
    int nPoints = 10;//(int) (200 + Ns / 200.0);     

    
    //#pragma omp parallel for
    for (int p = 0; p < nPartitions - 2; p++) {     

        if (p == 0) {
            std::cout << "begin parallel region " << omp_get_num_threads() << std::endl;  
        }   
        
        vector <double> range(sample.begin() + partitions[p], sample.begin() + partitions[p + 2] + 1); 

//        Block block = Block(range, (int) nPoints * range.size() / Ns, Ns, p);      
        Block block = Block(range, nPoints, Ns, p);

        #pragma omp critical
        {
            x.insert(x.end(), block.x.begin(), block.x.end()); 
        }    

        //#pragma omp critical
        blocks[p] = block;
        //blocks.push_back(block);
    }    
    sort(x.begin(), x.end());
}


void stitchPDF::stitch() {  
    
    double power = 2;
    
    double p0;
    double p1;
    double c0;
    double c1;
    double u0;
    double u1;
    
    int block0 = 0;
    int block1 = block0 + 1;
    double min = x[0];
    double max = blocks[0].xMax;
    double lowMax = blocks[block0].xMax;
    double cTotal;
    int nBlocks = blocks.size();
    if (nBlocks == 1) {
        x = blocks[0].x;
        pdf = blocks[0].pdf;
    } else {
        for (int i = 0; i < x.size(); i++) {         
            if (x[i] > lowMax) {
                if (++block0 > (nBlocks - 2)) {
                    block0--;
                    max = blocks[block1].xMax;
                } else {
                    max = blocks[block0].xMax;
                    block1 = block0 + 1;
                }
                lowMax = blocks[block0].xMax;
                min = blocks[block1].xMin;            
            }
            p0 = blocks[block0].pdfPoint(x[i]);
            p1 = blocks[block1].pdfPoint(x[i]);   
        
            u0 = (blocks[block0].cdfPoint(x[i]) - blocks[block0].cdfPoint(min)) /
                 (blocks[block0].cdfPoint(max) - blocks[block0].cdfPoint(min));
            u1 = (blocks[block1].cdfPoint(x[i]) - blocks[block1].cdfPoint(min)) /
                 (blocks[block1].cdfPoint(max) - blocks[block1].cdfPoint(min));
        
            c1 = pow(u1, power);
            if (x[i] > lowMax) {
                c0 = 0;
            } else {            
                c0 = pow(1 - u0, power);
            }        
            cTotal = c0 + c1;        
            pdf.push_back(p0 * c0 / cTotal + p1 * c1 / cTotal);       
        }
    }
//    WriteResults write;
//    write.writeColumns("stitched.txt", x, pdf, x.size());
}

void stitchPDF::branchNC() {
    int nBranch = 1;
    partitions.push_back(0);
    partitions.push_back(Ns - 1);
    for (int j = 1; j < maxLevel; j++) {
        bool addPartition = 0;
        for (int iBranch = 0; iBranch < nBranch; iBranch++) {
            vector <double> range(sample.begin() + partitions[iBranch], sample.begin() + partitions[iBranch + 1] + 1);
            int partition = uniformSplit(range);
            if (partition > 0) {
//                out.print("partition", partitions[iBranch] + partition);
                partitions.push_back(partitions[iBranch] + partition);
                addPartition = 1;
            }     
        }
        if (addPartition) {
            sort(partitions.begin(), partitions.end()); 
            nBranch = partitions.size() - 1;
        } else {
            break;
        }        
    }
}