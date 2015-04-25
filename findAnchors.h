#ifndef FINDANCHORS_H
#define FINDANCHORS_H
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <unordered_map>

#include "api/BamReader.h"
#include "api/BamWriter.h"

// for testing
#include <ctime>
#include <list>


typedef std::string barcode_str;
typedef std::string read_str;
typedef std::unordered_map <uint32_t, uint32_t> umapKmer;
typedef std::map<barcode_str, uint32_t> mapCount;           // counts the numbers of barcodes




// Question: any things special about handling pair end reads?
class readwKmer { 
    //FILE *_fa;       	     // filename frmo which the read comes from
    uint32_t _expLen; // length of the read
    uint32_t _skip;    // bases to skip to get the next kmer
    uint32_t _kmerLen; // length of the kmer

    uint32_t kmerpos;   // position of the kmer in the read
    read_str seq;      // sequence of the read
    int encodeKmer();
    int encode(char base);

public:
    barcode_str bc; // bar code of the read (BX tag in .bam file)
    bool eor; // whether the at the end of the read or not
    uint32_t readLen; // length of the read 
    bool clipped;   // if the read is clipped or not
    std::string kmer_str;    // kmer string to be converted to kmer
    uint32_t kmer;  // kmer in bits (k <= 16)
    uint32_t numAnchors;  // # of kmers anchored in the read
    std::map <uint32_t, uint32_t> matches; // pos in read -> pos in ref 
    bool anchored;  // true if at least one kmer in the read is anchored
    uint32_t anchorPos; // pos in ref where the read is anchored
    readwKmer(uint32_t expLen, uint32_t skip, uint32_t kmerLen); // constructor
    ~readwKmer(); // destructor

    int init(read_str str);
    int getNextKmer();
    int lookupKmer(umapKmer& m); 
    int printAll(void);

    	// look up the unique map, update matches, numAnchors, anchored
    int determineAnchor();
    	// use matches to determine the read is anchored and update anchorPos
};

// io.cpp
int file_to_unimap(const char* fname, umapKmer& m, const int k);
int file_to_wellmap(const char* fname, mapCount& m);
int map_to_file(const char* bFname, mapCount m, const int readsleft, const int areadsleft);


// printer.cpp - used for testing and displaying
int printbits(uint64_t binStr, int len);
int printbits32(uint32_t binStr);
int printbits64(uint64_t binStr);
int printMapinfo(umapKmer m);

#endif