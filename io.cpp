#include "findAnchors.h"


int file_to_unimap(const char* fname, umapKmer& m, const int k){
    std::ifstream file;
    uint64_t kmer;
    int64_t pos;
    uint64_t map_size;

    file.open(fname, std::ios::in | std::ios::binary);
    if (!file.is_open()) {
       printf("file_to_unimap: Cannot open the file %s!\n", fname);
       exit(1);
    }

    file.read(reinterpret_cast<char*>(&map_size), sizeof(map_size));
    for (uint64_t i = 0; i < map_size; i++) {
        file.read(reinterpret_cast<char*>(&kmer), sizeof(kmer));
        file.read(reinterpret_cast<char*>(&pos), sizeof(pos));
        if (pos > -1) {
            m[kmer] = (uint32_t) pos; // since <uint64_t, uint32_t> umapKmer 
        }
    }
    file.close();
    printf("Finished loading file %s!\n", fname);
    printf("Total number of distinct kmers:\t%lld\n", (long long int) map_size);
    printf("Distinct kmers that are unique:\t%lld\n", (long long int) m.size());
    return 0;
}

/*
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
*/