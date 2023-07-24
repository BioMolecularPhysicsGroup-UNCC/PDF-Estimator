/* 
 * File:   stitchPDF.cpp
 * Author: jfarm
 * 
 * Created on November 6, 2022, 10:25 AM
 */

#include "StitchPDF.h"

stitchPDF::stitchPDF(vector<double> sample, InputParameters input) {
    
    Ns = sample.size();    
    sort(sample.begin(), sample.end());       
    this->sample = sample;
    this->input = input;
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

    #pragma omp parallel for schedule(static)
    for (int b = 0; b < block_count; b++) {
        int start = b*block_size;
        int end = b == block_count - 1 ? sample.size() : (b+1)*block_size;

        branch(start, end - 1);
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
    double leftBoundary = left == 0 ? sample[0] : (sample[left] + sample[left - 1])/2;
    double rightBoundary = right == sample.size() - 1 
        ? sample[sample.size() - 1] 
        : (sample[right] + sample[right + 1])/2;

    // if leftBoundary is close enough (find a good epsilon as a hyperparameter) to rightBoundary,
    // we should compute the area the same (n / N), but how do we get the height? (pdf, derivative)
    // if we are in a degenerate block, find the first left-most and right-most non-too-close points and 
    // get their difference. this is dx

    if ((right - left >= 40) || !uniformSplit(left, right)) {

        LayerOptions layer = { true, left == 0, right == sample.size() - 1 };
        vector <double> range(sample.begin() + left, sample.begin() + right);
        Block block = Block(std::move(range), 10, Ns, left, out.debug, layer, input);

        // block.estimateBlock has side-effects and returns true on success -
        // short-circuiting is used to evaluate regardless of success and size
        if (block.estimateBlock(leftBoundary, rightBoundary) || (right - left <= 40)) {
            #pragma omp critical
            {
                partitions.push_back(left);
                blocks.push_back(block);
                //x.insert(x.end(), block.x.begin(), block.x.end());
            }
            return;
        }
    }
    
    branch(left, (right+left)/2);
    branch((right+left)/2 + 1, right);
}

bool stitchPDF::uniformSplit(int left, int right) {    
    int n = right - left;    

    if (n < 40) return false;

    double threshold = getThreshold(n);
    
    vector <double> dx;
    dx.reserve(n - 1);
    for (int i = left; i < right; i++) {
        dx.push_back(sample[i + 1] - sample[i]);
    }  
    
    double ratio = getRatio(dx);
    
    if (ratio < threshold) {
        return false;            
    } else {
out.print("ratio", ratio);
out.print("threshold", threshold);
    }
    return true;
}

inline double stitchPDF::getThreshold(int n) {
    static const double log10Nmin = log10((double)maxLevel);
    static const double log10Nmax = log10((double)maxBlock);

    double log10n = log10((double)n);
    double m = ( m0*(log10n - log10Nmin) + m5*(log10Nmax - log10n) )/(log10Nmax - log10Nmin);
    double b = ( b0*(log10n - log10Nmin) + b5*(log10Nmax - log10n) )/(log10Nmax - log10Nmin);
    return pow(10.0, m*log10n + b);
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


inline double stitchPDF::getRatio(const vector <double> &sample) {

    vector<double> maxes(sample.begin(), sample.begin() + window);
    vector<double> mins(sample.begin(), sample.begin() + window);

    re_min(mins, window);
    re_max(maxes, window);

    for (int j = window + 1; j < sample.size(); j++) {
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
            int left = ((partitions[l] + partitions[l + 1])/2) - BUFFER;
            int right = ((partitions[l + 1] + partitions[l + 2])/2) + BUFFER;

            vector<double> range(sample.begin() + left, sample.begin() + right); 
            LayerOptions layer = { false, false, false };

            blocks.insert(blocks.begin() + (2*l)+1, Block(std::move(range), 10, Ns, 2*(l+1), out.debug, layer, input));
        }

        #pragma omp parallel for schedule(static)
        for (int l = 0; l < partitions.size() - 2; l++) {
            // just recompute for now
            int left = ((partitions[l] + partitions[l + 1])/2) - BUFFER;
            int right = ((partitions[l + 1] + partitions[l + 2])/2) + BUFFER;

            double leftBoundary = (sample[left] + sample[left - 1])/2;
            double rightBoundary = (sample[right] + sample[right + 1])/2;

            blocks[(2*l)+1].estimateBlock(leftBoundary, rightBoundary);

            //#pragma omp critical
            //x.insert(x.end(), blocks[(2*l)+1].x.begin(), blocks[(2*l)+1].x.end());
        }
    }
}

void stitchPDF::stitch() {  

    //sort(x.begin(), x.end());
    
    double power = 2;
    
    double u1, u2, p1, p2, jMid, endMid, jMidInK, jMaxInK, jMinInK, jMinCdf;

    int j = 0; // bottom block
    int k = 1; // top block
    int ji = 0;
    int ki = 0;
    
    int nBlocks = blocks.size();

    int xSize = 0;

    for (int i = 0; i < blocks.size(); i++) {
        xSize += blocks[i].x.size();
    }

    if (nBlocks == 1) {
        x = blocks[0].x;
        pdf = blocks[0].pdf;
        cdf = blocks[0].cdf;
    } else {
#ifdef _OLD_PDF_STITCH
        // get the last two x values and store them for computing CDF.
        double prev = blocks[0].x[0];
        double earlierPrev = 0;

        double topBlockStart = blocks[1].xMin;
        endMid = (blocks[blocks.size() - 1].xMin + blocks[blocks.size() - 1].xMax) / 2;

        jMid = (blocks[j].xMin + blocks[j].xMax) / 2;
        jMidInK = blocks[k].cdfPoint(jMid);
        jMaxInK = blocks[k].cdfPoint(blocks[j].xMax);
        jMinInK = blocks[k].cdfPoint(blocks[j].xMin);
        jMinCdf = blocks[j].cdfPoint(jMid);

        pdf.resize(xSize);

        bool indicator = false;

        //#pragma omp parallel for schedule(static) firstprivate(j, k, indicator, jMid)
        for (unsigned int i = 0; i < xSize; i++) {
            /*
            if (!indicator) {
                indicator = true;

                while (x[i] > blocks[j].xMax) j += 2;

                k = (j == blocks.size() - 1 || (j > 0 && x[i] < blocks[j - 1].xMax)) 
                    ? j - 1 
                    : j + 1;

                jMid = (blocks[j].xMin + blocks[j].xMax) / 2;
            }
            */

            if (prev > endMid) {
                earlierPrev = prev;
                prev = pdf[i] = blocks[blocks.size() - 1].x[ji++];
            } else if (prev < topBlockStart) {
                earlierPrev = prev;
                prev = pdf[i] = blocks[0].x[ji++];  
            } else {
                if (j < k) {
                    u1 = (1 - blocks[j].cdfPointHint(prev, ji)) / 
                        (1 - jMidInK);

                    u2 = (blocks[k].cdfPointHint(prev, ki) - jMidInK) /
                        (jMaxInK - jMidInK);

                    p1 = blocks[j].pdfPointHint(prev, ji) * u1;
                    p2 = blocks[k].pdfPointHint(prev, ki) * u2;
                } else {
                    u1 = ((1 - blocks[k].cdfPointHint(prev, ki)) - (1 - jMidInK)) / 
                        ((1 - jMinInK) - (1 - jMidInK));

                    u2 = blocks[j].cdfPoint(prev) / jMinCdf;

                    p1 = blocks[k].pdfPointHint(prev, ki) * u1;
                    p2 = blocks[j].pdfPointHint(prev, ji) * u2;
                }   

                pdf[i] = (p1 + p2) / (u1 + u2);

                double nextJ = ji < blocks[j].x.size() - 1 ? blocks[j].x[ji + 1] : blocks[j + 2].x[0];
                double nextK = k < blocks.size() - 2 
                    ? (blocks[k].x[ki + 1] < jMid ? blocks[k].x[ki + 1] : blocks[k + 2].x[BUFFER])
                    : blocks[k].x[ki];


                earlierPrev = prev;
                if (nextJ < nextK) {
                    prev = nextJ;
                    if (ji == blocks[j].x.size() - 1) {
                        ji = 0;
                        j += 2;

                        jMid = (blocks[j].xMin + blocks[j].xMax) / 2;
                        jMinCdf = blocks[j].cdfPoint(jMid);
                        jMidInK = blocks[k].cdfPoint(jMid);
                        jMaxInK = blocks[k].cdfPoint(blocks[j].xMax);
                        jMinInK = blocks[k].cdfPoint(blocks[j].xMin);
                    } else {
                        ji++;
                    }
                } else {
                    prev = nextK;
                    if (ki == blocks[k].x.size() - 1) {
                        ki = BUFFER;
                        k += 2;

                        jMidInK = blocks[k].cdfPoint(jMid);
                        jMaxInK = blocks[k].cdfPoint(blocks[j].xMax);
                        jMinInK = blocks[k].cdfPoint(blocks[j].xMin);
                    } else {
                        ki++;
                    }
                }
            }

            if (i != 0) {
                double average = 0.5 * (pdf[i] + pdf[i-1]);
                double area = average * (prev - earlierPrev);
                //cdf[i] = cdf[i - 1] + area;  
            }
        }
#else

        int arrsize = (blocks.size() + 1) * input.resolution;
        x.resize(arrsize);
        pdf.resize(arrsize);
        cdf.resize(arrsize);

        double start, end, step,
                bottomExtremaInTop, bottomMidInTop, bottomMidCdf;

        #pragma omp parallel for schedule(static)
        for (int i = 0; i <= blocks.size(); i++) {
            int last;

            if (i == blocks.size()) {
                start = (blocks[i - 1].xMin + blocks[i - 1].xMax) / 2;
                end = blocks[i - 1].xMax;
            } else {
                if (i % 2 == 0) {
                    start = blocks[i].xMin;
                    end = (blocks[i].xMin + blocks[i].xMax) / 2;
                    
                    if (i > 0) {
                        bottomMidInTop = blocks[i - 1].cdfPoint(end);
                        bottomExtremaInTop = blocks[i - 1].cdfPoint(start);
                    }
                } else {
                    start = (blocks[i - 1].xMin + blocks[i - 1].xMax) / 2;
                    end = blocks[i - 1].xMax;

                    bottomMidInTop = blocks[i].cdfPoint(start);
                    bottomExtremaInTop = blocks[i].cdfPoint(end);
                    bottomMidCdf = blocks[i - 1].cdfPoint(start);
                }
            }

            step = (end - start)/((double)input.resolution + 1.0);

            for (int j = 0; j < input.resolution; j++) {

                double xPoint = start + (step * (j + 1));
                x[(input.resolution * i) + j] = xPoint;

                // loop unswitching should pull this out to prevent repeated conditional checks
                if (i == 0 || i == blocks.size()) {
                    pdf[(input.resolution * i) + j] = blocks[min((size_t)i, blocks.size() - 1)].pdfPoint(xPoint);
                } else if (i % 2 == 1) {
                    u1 = (1 - blocks[i - 1].cdfPoint(xPoint)) / 
                        (1 - bottomMidInTop);

                    u2 = (blocks[i].cdfPoint(xPoint) - bottomMidInTop) /
                        (bottomExtremaInTop - bottomMidInTop);

                    p1 = blocks[i - 1].pdfPoint(xPoint) * u1;
                    p2 = blocks[i].pdfPoint(xPoint) * u2;

                    pdf[(input.resolution * i) + j] = (p1 + p2) / (u1 + u2);
                } else {
                    u1 = ((1 - blocks[i - 1].cdfPoint(xPoint)) - (1 - bottomMidInTop)) / 
                        ((1 - bottomExtremaInTop) - (1 - bottomMidInTop));

                    u2 = blocks[i].cdfPoint(xPoint) / bottomMidCdf;

                    p1 = blocks[i - 1].pdfPoint(xPoint) * u1;
                    p2 = blocks[i].pdfPoint(xPoint) * u2;

                    pdf[(input.resolution * i) + j] = (p1 + p2) / (u1 + u2);
                }
            }
        }
#endif

        cdf[0] = 0.0;

        for (int i = 1; i < pdf.size(); i++) {
            cdf[i] = cdf[i - 1] + ((0.5 * (pdf[i] + pdf[i - 1])) * (x[i] - x[i - 1]));
        }
    }
    
    
    
//    WriteResults write;
//    write.writeColumns("stitched.txt", x, pdf, x.size());
}
