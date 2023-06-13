

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
    out.print("branch");
    branch(0, Ns - 1);  
//    branchNC();
    out.print("stagger");
    stagger();
    out.print("estiimate"); 
    estimate(); 
    out.print("stitch");  
    stitch(); 
}

stitchPDF::stitchPDF(const stitchPDF& orig) {
}

stitchPDF::~stitchPDF() {
}

vector <double> stitchPDF::getX(int nPoints){
    WriteResults write;
    vector <double> gridPoints;
    double dx = 1.0 / (nPoints - 1);
    double dxSum = 0;
    for (int i = 0; i < nPoints; i++) {
        gridPoints.push_back(dxSum);
        dxSum += dx;
    }        
    return write.interpolateGrid(cdf, x, gridPoints);    
}

vector <double> stitchPDF::getPDF(int nPoints){    
    WriteResults write;
    vector <double> gridPoints;
    double dx = 1.0 / (nPoints - 1);
    double dxSum = 0;
    for (int i = 0; i < nPoints; i++) {
        gridPoints.push_back(dxSum);
        dxSum += dx;
    }        
    return write.interpolateGrid(cdf, pdf, gridPoints);   
    
}

vector <double> stitchPDF::getCDF(int nPoints){
    vector <double> gridPoints;
    double dx = 1.0 / (nPoints - 1);
    double dxSum = 0;
    for (int i = 0; i < nPoints; i++) {
        gridPoints.push_back(dxSum);
        dxSum += dx;
    }        
    return gridPoints;     
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
        out.print("partition ", partition);
        if ((right - partition) > window) {
            branch(partition, right);
        }
    }
}


int stitchPDF::uniformSplit(const vector <double> &sample) {    
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
    } else {
out.print("ratio", average);
out.print("threshold", threshold);
    }
    return partition;
}

static inline void re_max(vector<double> &maxes, int window) {
    double curr_least_max = maxes[0];
    int curr_least_max_index = 0;

    for (int j = 1; j < window; j++) {
        if (maxes[j] < curr_least_max) {
            curr_least_max = maxes[j];
            curr_least_max_index = j;
        }
    }

    std::swap(maxes[0], maxes[curr_least_max_index]);
}

static inline void re_min(vector<double> &mins, int window) {
    double curr_greatest_min = mins[0];
    int curr_greatest_min_index = 0;

    for (int j = 1; j < window; j++) {
        if (mins[j] > curr_greatest_min) {
            curr_greatest_min = mins[j];
            curr_greatest_min_index = j;
        }
    }

    std::swap(mins[0], mins[curr_greatest_min_index]);
}


double stitchPDF::getRatio(const vector <double> &sample, int start, int stop) {
    int n = stop - start + 1;    

    vector<double> maxes(sample.begin() + start, sample.begin() + start + window);
    vector<double> mins(sample.begin() + start, sample.begin() + start + window);

    re_min(mins, window);
    re_max(maxes, window);

    for (int j = start + window + 1; j < stop; j++) {
        if (sample[j] > maxes[0]) {
            maxes[0] = sample[j];
            re_max(maxes, window);
        } 
        
        if (sample[j] < mins[0]) {
            mins[0] = sample[j];
            re_min(mins, window);
        }
    }

    double dxMin = 0;
    double dxMax = 0;
    for (int i = 0; i < window; i++) {
        dxMin += mins[i];
        dxMax += maxes[i];
    }       

    return dxMin == 0 ? 0 : dxMax/dxMin;    
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
    int nPoints = 10;//(int) (200 + Ns / 200.0);      
    
//    parallel_for(nPartitions - 2, [&](int start, int end){ 
//        for(int p = start; p < end; p++) {            
    for (int p = 0; p < nPartitions - 2; p++) {            
        vector <double> range(sample.begin() + partitions[p], sample.begin() + partitions[p + 2] + 1);        
//        Block block = Block(range, (int) nPoints * range.size() / Ns, Ns, p);      
        Block block = Block(range, nPoints, Ns, p, out.debug);
        x.insert(x.end(), block.x.begin(), block.x.end());
        blocks.push_back(block);       
    }            
//    }
//);
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
        for (unsigned int i = 0; i < x.size(); i++) {         
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
    cdf.push_back(0);
    for (unsigned int i = 1; i < x.size(); i++) {   
        double average = 0.5 * (pdf[i] + pdf[i-1]);
        double area = average * (x[i] - x[i - 1]);
        cdf.push_back(cdf[i - 1] + area);  
    }
    
    
//    WriteResults write;
//    write.writeColumns("stitched.txt", x, pdf, x.size());
}

void stitchPDF::branchNC() {
    int nBranch = 1;
    partitions.push_back(0);
    partitions.push_back(Ns - 1);
    for (int j = 1; j < maxLevel; j++) {
        out.print("j", j);
        bool addPartition = 0;
        for (int iBranch = 0; iBranch < nBranch; iBranch++) {
            vector <double> range(sample.begin() + partitions[iBranch], sample.begin() + partitions[iBranch + 1] + 1);
            int partition = uniformSplit(range);
            if (partition > 0) {
//                out.print("partition", partitions[iBranch] + partition);
                partitions.push_back(partitions[iBranch] + partition);
                out.print("partition", partitions[iBranch] + partition);
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
