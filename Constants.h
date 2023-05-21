#ifndef _CONSTANTS_PDFESTIMATOR
#define _CONSTANTS_PDFESTIMATOR 1

// If the input is bigger than this, split the file up among MPI processes and compute in parallel
#define MAX_SERIAL_FILE_SIZE 300

// No block should be more than 50,000 * 2 = 100,000 elements wide
#define MAX_2_BLOCK_SIZE 50000

// If our file has more than this, we need to revert to chunking mechanism
// Also, if the size of the file divided by the number of processes is bigger than this,
// we cannot reliably fit the data into memory, so fail.
#define MAX_CHUNK_SAMPLE_SIZE (4294967296/8)


#endif