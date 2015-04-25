#include "findAnchors.h"

readwKmer::readwKmer(uint32_t expLen, uint32_t skip, uint32_t kmerLen) {
	// constructor
	// set up the private parameters
	_expLen = expLen;
	_skip = skip;
	_kmerLen = kmerLen;
	// initialize public parameters
	// initialize map?
}

readwKmer::~readwKmer(void) {
	// destructor
	// printf("Object is deleted.\n");
}

int readwKmer::init(read_str str) {
	kmerpos = 0;
	numAnchors = 0;
	anchored = false;
	anchorPos = 0; 
	matches.clear(); 
	seq.assign(str);
	readLen = str.length();
	// check is the read is hard clipped 
	if (readLen < _expLen) {
		clipped = true;
	} else {
		clipped = false;
	}
	// catch the case when the read is clipped too much
	if (readLen < _kmerLen) {
		eor = true;
	} else {
		eor = false;
		encodeKmer();
	}
	return 0;
}

int readwKmer::encodeKmer(void) {
	// TO DO: speed up this part of the code if necessary
	// use kmerpos and seq to get the kmer in bits
	kmer = 0; // initialize the kmer to be zero
	kmer_str = seq.substr(kmerpos,_kmerLen);
	for (uint32_t i = 0; i < _kmerLen; ++i) { 
		kmer =  kmer << 2 | encode(kmer_str.at(i));
	}
	return 0;
}

int readwKmer::encode(char base) {
    if      (base == 'T' || base == 't') {return 3;} // 11
    else if (base == 'G' || base == 'g') {return 2;} // 10
    else if (base == 'C' || base == 'c') {return 1;} // 01
    else {return 0;} // A and other letters chars are all 0 (but we have hasN to control)
}

int readwKmer::getNextKmer(void) {
	// skip the number of bases to get the next kmer if eor is false
	kmerpos = kmerpos + _skip;
	if (kmerpos > readLen - _kmerLen) { // kmerpos near the end
		eor = true;
	} 
	else {
		encodeKmer();
	}
	return 0;
}

int readwKmer::lookupKmer(umapKmer& m) {
	if (m.count(kmer) > 0) { // if the kmer is a key in the unique kmer map
		anchored = true; // found an anchor
		numAnchors ++;
		matches[kmerpos] = m[kmer]; // kmer pos in read -> position in reference
	}
	return 0;
}


int readwKmer::determineAnchor(void) {
	// TO DO: handle matches to update anchorPos 

	return 0;
}

int readwKmer::printAll(void) {
	//printf("%u\t", kmerpos);
	//printbits32(kmer);
	printf("\t......%s\t", kmer_str.c_str());
	printf("\t number of anchors: %u\n", numAnchors);
	
	// print the map of matches
	uint32_t kpos;
    uint32_t val;
    std::map<uint32_t, uint32_t>::iterator it;
    for (it = matches.begin(); it != matches.end(); ++it) {
        kpos = it->first;
        val = it->second;
        printf("\tread position: %u,", kpos);
        printf(" anchored at reference position: %u\n", val);
    }
	return 0;
}