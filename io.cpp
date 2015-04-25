#include "findAnchors.h"


int file_to_unimap(const char* fname, umapKmer& m, const int k){
	std::ifstream file;
    uint64_t kmer64;
    uint32_t kmer;
    uint64_t pos64;
    uint64_t pos;
    uint64_t map_size;


	file.open(fname, std::ios::in | std::ios::binary);
	if (!file.is_open()) {
       printf("file_to_unimap: Cannot open the file %s!\n", fname);
       exit(1);
    }

    if (k > 16){
        printf("file_to_unimap: k must be no more than 16");
        exit(1);
    }

    file.read(reinterpret_cast<char*>(&map_size), sizeof(map_size));
    for (uint64_t i = 0; i < map_size; i++) {
    	file.read(reinterpret_cast<char*>(&kmer64), sizeof(kmer64));
        file.read(reinterpret_cast<char*>(&pos64), sizeof(pos64));
        kmer = kmer64 & 0xFFFFFFFF; // convert 64 to 32
        pos = pos64 & 0xFFFFFFFF; // convert 64 to 32
        m[kmer] = pos;
    }
	file.close();
    printf("Finished loading file %s!\n", fname);
    return 0;
}

int file_to_wellmap(const char* fname, mapCount& m) {
    std::string line;
    std::ifstream file (fname);
    if (!file.is_open()) {
       printf("file_to_map: Cannot open the file %s!\n", fname);
       exit(1);
    }

    while (std::getline (file,line)) {
        m[line] = 0;
    }
    file.close();
    printf("Finished loading file %s!\n", fname);
    return 0;
}

int map_to_file(const char* bFname, mapCount m, const int readsleft, const int areadsleft) {
    std::string fname1(bFname);   
    fname1 += std::string("_wacount"); 
    std::ofstream file1 (fname1.c_str());
    if (!file1.is_open()) {
       printf("map_to_file: Cannot open the file %s!\n", fname1.c_str());
       exit(1);
    }
    std::string fname3(bFname);   
    fname3 += std::string("_wasummary"); 
    std::ofstream file3 (fname3.c_str());
    if (!file3.is_open()) {
       printf("map_to_file: Cannot open the file %s!\n", fname3.c_str());
       exit(1);
    }

    uint32_t wellswanchors = 0;
    std::string bc; // barcode
    uint32_t ct;    // count
    for ( mapCount::iterator it = m.begin(); it != m.end(); ++it) {
        bc = it->first;
        ct = it->second;
        file1 << ct << "\n"; 

        if (ct > 0) {
            wellswanchors ++;
        }
    }

    file3 << "number of wells with anchors: " <<  wellswanchors << " / " << m.size() <<"\n";
    file3 << "number of leftover reads anchored: " << areadsleft << "/ " << readsleft << "\n";

    return 0;
}