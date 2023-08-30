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
#include <sstream>
#include <limits>
#include "InputParameters.h"
#include "Variable.h"
#include "JointProbability.h"
#include "InputData.h"
#include "ScoreQZ.h"
#include "MinimizeScore.h"
#include "WriteResults.h"
#include "StitchPDF.h"
#include "OutputControl.h"

using namespace std;

#ifdef outputCommandLine

int main(int argc, char** argv) {               
    bool stitch = true;
    double readTime;                                                  
    readTime = time(NULL);
    OutputControl out;   
    out.debug = true;

   
    #ifndef PDF_BINARY_1D_READ
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
            InputParameters input;
            stitchPDF stitch = stitchPDF(samples, input);                          
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
    #else
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    InputParameters input;   
    if (!input.userInput(argc, argv)) {
        return false;
    }

    MPI_Offset file_size;
    MPI_File sample_file;
    MPI_Status read_status;
    vector<double> samples;
    int* receive_sizes = new int[size];

    MPI_File_open(MPI_COMM_WORLD, (input.inputPath + input.inputFile).c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &sample_file);
    MPI_File_get_size(sample_file, &file_size);

    char* file_data;

    if (file_size > MAX_SERIAL_FILE_SIZE) {

        file_data = (char*)malloc(sizeof(char) * ((file_size/size) + 100));

        MPI_File_read_at_all(sample_file, (file_size/size)*rank, file_data, (file_size/size) + 100, MPI_CHAR, &read_status);

        string file_string(file_data, (file_size/size) + 100);
        istringstream file_stream(file_string);

        string line;

        // discard potentially truncated entry
        if (rank != 0) {
            getline(file_stream, line);
        }

        while (file_stream.tellg() < (file_size/size) + 1 && file_stream.tellg() != -1) {
            getline(file_stream, line);
            double sample = stod(line);
            samples.push_back(sample);
        }

        free(file_data);

        int sample_size = samples.size();
        MPI_Offset total_sample_size = 0;

        MPI_Allgather(&sample_size, 1, MPI_INT, receive_sizes, 1, MPI_INT, MPI_COMM_WORLD);

        for (int i = 0; i < size; i++) {
            total_sample_size += receive_sizes[i];
        }

        sort(samples.begin(), samples.end());
        vector<double> new_samples;

        int division = 1;
        int* rank_sizes = new int[size];

        std::copy(receive_sizes, receive_sizes + size, rank_sizes);

        bool sent = false;
        bool oversized = false;
    
        while ((division *= 2, division) < size*2) {
            int* temp_ranks = new int[size];
            std::copy(rank_sizes, rank_sizes + size, temp_ranks);

            for (int i = 0; i < size - division; i += division) {
                long int sum = (long int)rank_sizes[i] + (long int)rank_sizes[i + division/2];

                if (sum > MAX_CHUNK_SAMPLE_SIZE) {
                    oversized = true;
                    std::copy(temp_ranks, temp_ranks + size, rank_sizes);
                    break;
                } else {
                    rank_sizes[i] = (int)sum;
                }
            }

            if (oversized) {
                break;
            } else if (sent) {
                continue;
            }

            if (rank % division == 0) { // We are receiving

                if (rank + (division/2) >= size) {
                    continue;
                }

                int receive_size = rank_sizes[rank + (division/2)];

                samples.resize(sample_size + receive_size);
                MPI_Recv(&samples[sample_size], receive_size, MPI_DOUBLE, rank + (division/2), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                vector<double> new_samples;

                int left_index = 0;
                int right_index = sample_size;

                for (int i = 0; i < sample_size + receive_size; i++) {
                    if ((samples[left_index] < samples[right_index] && left_index < sample_size) || right_index == sample_size + receive_size) {
                        new_samples.push_back(samples[left_index]);
                        left_index++;
                    } else {
                        new_samples.push_back(samples[right_index]);
                        right_index++;
                    }
                }

                // this should fully dealloc the old samples vector
                samples = new_samples;
                sample_size += receive_size;

            } else { // We are sending
                sent = true;
                MPI_Send(&samples[0], sample_size, MPI_DOUBLE, rank - (division/2), 0, MPI_COMM_WORLD);
            }
        }

        if (!oversized && rank == 0) {
            ofstream fileout("temp.txt", ios::out);

            for (double sample : samples) {
                fileout << std::fixed << sample << std::endl;
            }

            fileout.close();
        } else if (oversized) {
            division /= 2;

            // we couldn't fit everything into a single node, we need to chunk
            int chunk_size = (MAX_CHUNK_SAMPLE_SIZE - rank_sizes[0])/(size/division);

            int divided_size = division == 1 ? size : (size/division) + 1;

            int* active_sizes = new int[divided_size];
            int* incoming_sizes = new int[divided_size];
            int* displacements = new int[divided_size];
            int* offset = new int[divided_size];

            displacements[0] = 0;
            incoming_sizes[0] = rank_sizes[0];

            for (int k = 0; k < size; k += division) {
                active_sizes[k/division] = rank_sizes[k];
                offset[k] = 0;

                if (k != 0) {
                    incoming_sizes[k/division] = std::min(chunk_size, rank_sizes[k]);
                    displacements[k/division] = displacements[(k-division)/division] + incoming_sizes[(k-division)/division];
                }
            }

            if (rank == 0) {
                samples.resize(MAX_CHUNK_SAMPLE_SIZE);

                MPI_Gatherv(MPI_IN_PLACE, sample_size, MPI_DOUBLE, &samples[0], 
                    incoming_sizes, displacements, MPI_DOUBLE, 0, MPI_COMM_WORLD);

                for (int k = division; k < size; k += division) {               
                    active_sizes[k/division] -= incoming_sizes[k/division];
                }

                vector<double> new_samples;

                // to be configured, just temporary for now
                ofstream binout("temp.bin", ios::out | ios::binary);

                for (MPI_Offset i = 0; i < total_sample_size; i++) {
                    double minimum = numeric_limits<double>::max();
                    int index = 0;

                    for (int j = 0; j < divided_size; j++) {
                        if (incoming_sizes[j] != 0 && samples[offset[j] + displacements[j]] < minimum) {
                            minimum = samples[offset[j] + displacements[j]];
                            index = j;
                        }
                    }

                    new_samples.push_back(minimum);
                    offset[index]++;

                    if (offset[index] == incoming_sizes[index]) {

                        if (index) {
                            incoming_sizes[index] = std::min(chunk_size, active_sizes[index]);
                            active_sizes[index] -= incoming_sizes[index];

                            if (incoming_sizes[index] != 0) 
                                MPI_Recv(&samples[displacements[index]], incoming_sizes[index], MPI_DOUBLE, index * division, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                            
                            offset[index] = 0;
                        } else {
                            incoming_sizes[0] = 0;
                        }
                    }

                    if (new_samples.size() == MAX_CHUNK_SAMPLE_SIZE) {
                        binout.write((char*)&new_samples[0], sizeof(double) * MAX_CHUNK_SAMPLE_SIZE);
                        new_samples.clear();
                        // should chunk here
                    }

                }

                binout.write((char*)&new_samples[0], sizeof(double) * new_samples.size());
                new_samples.clear();
                binout.close();
                
            } else if (rank % division == 0) {

                MPI_Gatherv(&samples[0], incoming_sizes[rank/division], MPI_DOUBLE, NULL, 
                        NULL, NULL, MPI_DOUBLE, 0, MPI_COMM_WORLD);

                int remaining = sample_size - incoming_sizes[rank/division];

                while (remaining > 0) {
                    MPI_Ssend(&samples[sample_size - remaining], std::min(chunk_size, remaining), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
                    remaining -= std::min(chunk_size, remaining);
                }
            }
        } 

        std::cout << samples.size() << " compared with " << total_sample_size << std::endl;
    } else {
        if (rank == 0) {
            file_data = (char*)malloc(sizeof(char) * file_size);

            MPI_File_read(sample_file, file_data, file_size, MPI_CHAR, MPI_STATUS_IGNORE);

            string file_string(file_data, file_size);
            istringstream file_stream(file_string);

            string line;

            while (getline(file_stream, line)) {
                double sample = stod(line);
                samples.push_back(sample);
            }

            free(file_data);
            file_stream.clear();
            file_string.clear();

            sort(samples.begin(), samples.end());

            //ofstream fileout("temp2.bin", ios::out | ios::binary);
            //fileout.write((char*)&samples[0], samples.size() * sizeof(double));
            //fileout.close();

            ofstream fileout("temp2.txt", ios::out);

            for (double sample : samples) {
                fileout << std::fixed << sample << std::endl;
            }

            fileout.close();

            std::cout << samples.size() << " in serial" << std::endl;
        }
    }

    MPI_File_close(&sample_file);

    MPI_Finalize();

    #endif
    return 0;
}
#endif