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
    /*
    out.print("branch");
    branch(0, Ns - 1);  
//    branchNC();
    out.print("stagger");
    stagger();
    out.print("estiimate"); 
    estimate(); 
    out.print("stitch");  
    stitch(); 
    */
    size_t block_count = std::max((size_t)1, sample.size()/50000);
    int block_size = sample.size()/block_count;

    #pragma omp parallel for
    for (int b = 0; b < block_count; b++) {
        int start = b*block_size;
        int end = b == block_count - 1 ? sample.size() - 1 : (b+1)*block_size;

        branch(start, end);
    }

    // Sort the blocks by their placement 
    std::sort(blocks.begin(), blocks.end(), [](Block lhs, Block rhs) { return lhs.blockNumber < rhs.blockNumber; });
    std::sort(partitions.begin(), partitions.end());
    partitions.push_back(sample.size());

    // Layer one blocks are numbered 2k + 1
    for (int i = 0; i < blocks.size(); i++) {
        blocks[i].blockNumber = 2*i + 1;
    }

    stagger();
    
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
    #pragma omp parallel 
    {
        vector <double> range(sample.begin() + left, sample.begin() + right);

        double leftBoundary = left == 0 ? sample[0] : (sample[left] + sample[left - 1])/2;
        double rightBoundary = right == sample.size() - 1 
            ? sample[sample.size() - 1] 
            : (sample[right] + sample[right + 1])/2;

        Block block = Block(range, 10, Ns, left, out.debug, true);

        if (uniformSplit(left, right) || (!block.estimateBlock(leftBoundary, rightBoundary) && (right - left >= 40))) {
            #pragma omp task untied final(right - left < 1000)
            branch(left, (right+left)/2);

            branch((right+left)/2 + 1, right);

        } else {
            #pragma omp critical
            {
                partitions.push_back(left);
                blocks.push_back(block);
                x.insert(x.end(), block.x.begin(), block.x.end());
            }
        }

        #pragma omp taskyield
    }
}


bool stitchPDF::uniformSplit(int left, int right) {    
    int n = right - left;    

    if (n < 40) return false;

    double threshold = thresholdCoeff * pow(n, thresholdExp);    
    int partition = ceil(n/2);
    
    vector <double> dx;
    dx.reserve(n - 1);
    for (int i = left; i < right; i++) {
        dx.push_back(sample[i + 1] - sample[i]);
    }  
    
    double ratioLeft = getRatio(dx, 0, partition - 1);
    double ratioRight = getRatio(dx, partition - 1, n - 1);     
    double average = (ratioLeft + ratioRight) / 2;
    if (average < threshold) {
        return false;            
    } else {
out.print("ratio", average);
out.print("threshold", threshold);
    }
    return true;
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
    if (blocks.size() > 1) {
        for (int l = 0; l < partitions.size() - 2; l++) {
            int left = ((partitions[l] + partitions[l + 1])/2) - 8;
            int right = ((partitions[l + 1] + partitions[l + 2])/2) + 8;

            vector<double> range(sample.begin() + left, sample.begin() + right); 

            blocks.insert(blocks.begin() + (2*l)+1, Block(range, 10, Ns, 2*(l+1), out.debug, false));
        }

        #pragma omp parallel for
        for (int l = 0; l < partitions.size() - 2; l++) {
            // just recompute for now
            int left = ((partitions[l] + partitions[l + 1])/2) - 8;
            int right = ((partitions[l + 1] + partitions[l + 2])/2) + 8;

            double leftBoundary = (sample[left] + sample[left - 1])/2;
            double rightBoundary = (sample[right] + sample[right + 1])/2;

            blocks[(2*l)+1].estimateBlock(leftBoundary, rightBoundary);

            #pragma omp critical
            x.insert(x.end(), blocks[(2*l)+1].x.begin(), blocks[(2*l)+1].x.end());
        }
    }
}

void stitchPDF::stitch() {  

    sort(x.begin(), x.end());
    
    double power = 2;
    
    double u1, u2, p1, p2, jMid, endMid;

    int j = 0; // bottom block
    int k = 1; // top block
    
    int nBlocks = blocks.size();

    if (nBlocks == 1) {
        x = blocks[0].x;
        pdf = blocks[0].pdf;
    } else {
        double topBlockStart = blocks[1].xMin;
        endMid = (blocks[blocks.size() - 1].xMin + blocks[blocks.size() - 1].xMax) / 2;
        jMid = (blocks[j].xMin + blocks[j].xMax) / 2;

        pdf.resize(x.size());

        bool indicator = false;

        #pragma omp parallel for schedule(static) firstprivate(j, k, indicator, jMid)
        for (unsigned int i = 0; i < x.size() - 1; i++) {
            if (!indicator) {
                indicator = true;

                while (x[i] < blocks[j].xMin) j += 2;

                k = (j > 0 && x[i] < blocks[j - 1].xMax) ? j - 1 : j + 1;

                jMid = (blocks[j].xMin + blocks[j].xMax) / 2;
            }

            if (x[i] > endMid)
                pdf[i] = x[i];
            else if (x[i] < topBlockStart) {
                pdf[i] = x[i];  
            } else {
                if (j < k) {
                    u1 = (1 - blocks[j].cdfPoint(x[i])) / 
                        (1 - blocks[j].cdfPoint(jMid));

                    u2 = (blocks[k].cdfPoint(x[i]) - blocks[k].cdfPoint(jMid) ) /
                        (blocks[k].cdfPoint(blocks[j].xMax) - blocks[k].cdfPoint(jMid));

                    p1 = blocks[j].pdfPoint(x[i]) * u1;
                    p2 = blocks[k].pdfPoint(x[i]) * u2;

                    if (x[i + 1] > blocks[j].xMax) {
                        j += 2;

                        jMid = (blocks[j].xMin + blocks[j].xMax) / 2;
                    }

                } else {
                    u1 = ((1 - blocks[k].cdfPoint(x[i])) - (1 - blocks[k].cdfPoint(jMid))) / 
                        ((1 - blocks[k].cdfPoint(blocks[j].xMin)) - (1 - blocks[k].cdfPoint(jMid)));

                    u2 = blocks[j].cdfPoint(x[i]) / blocks[j].cdfPoint(jMid);

                    p1 = blocks[k].pdfPoint(x[i]) * u1;
                    p2 = blocks[j].pdfPoint(x[i]) * u2;

                    if (x[i + 1] > jMid && x[i] < jMid) {
                        k += 2;
                    }
                }   

                pdf[i] = (p1 + p2) / (u1 + u2);
            }
        }
    }

    pdf[x.size() - 1] = x[x.size() - 1];

    cdf.push_back(0);
    for (unsigned int i = 1; i < x.size(); i++) {   
        double average = 0.5 * (pdf[i] + pdf[i-1]);
        double area = average * (x[i] - x[i - 1]);
        cdf.push_back(cdf[i - 1] + area);  
    }
    
    
//    WriteResults write;
//    write.writeColumns("stitched.txt", x, pdf, x.size());
}
