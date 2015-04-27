#include "findAnchors.h"

int main(int argc, char* argv[]){
    char* bamFname;
    char* mapFname;
    int k;          // kmer length
    int rLen;       // read length
    int intvl;      // # of bases to skip in read
    char* tName;    // BX

    if ( argc == 7 ) 
    {
        bamFname  = argv[1];
        mapFname = argv[2];
        k = atoi(argv[3]); 
        rLen = atoi(argv[4]);
        intvl = atoi(argv[5]); 
        tName = argv[6]; 
    } 
    else 
    {
        std::cerr << "Wrong number of arguments." << std::endl;
        exit(1);
    }

    // -------------- load the unique kmers into an unordered map --------------
    umapKmer uniqueKmers; // unique kmer -> string    
    file_to_unimap(mapFname, uniqueKmers, k); // load uniqueKmerMap (unordered)

    
    // -------------- FILE I/O --------------
    std::string bname(bamFname);
    std::string ofname(bamFname);
    ofname += std::string("_anc"); 
    ofname += std::to_string(static_cast<long long>(k)); 
    
    std::ofstream file (ofname.c_str());
    if (!file.is_open()) 
    {
        printf("Cannot open the file %s!\n", ofname.c_str());
        exit(1);
    }
    BamTools::BamReader reader;
    if (!reader.Open(bname)) 
    {
        printf("Cannot open bam file %s!\n", bname.c_str());
        exit(1);
    }


    // read the bam file and find anchored reads
    BamTools::BamAlignment al;
    readwKmer read(rLen, intvl, k);
    uint32_t totalR = 0;        // total # of reads
    uint32_t totalAnchored = 0; // total # of anchored reads
    uint32_t anchoredR = 0;     // # of anchored reads per barcode
    uint32_t loReads = 0;       // # of reads w/o barcodes
    uint32_t loAnchored = 0;    // # of anchored reads w/o barcodes
    uint32_t clipcount = 0;     // # of reads that are hard clipped
    barcode_str bc, bc_prev;
    char tagTypeName;           // type of tag field in BAM format
    bool nobc_flag = true;
    double duration;
    std::clock_t start;


    while (reader.GetNextAlignment(al))   // each BAM entry is processed in this loop
    { 
        read.init(al.QueryBases);
            
        while (!read.eor) 
        {   
            read.lookupKmer(uniqueKmers);
            read.getNextKmer();
        }

        read.determineAnchor(); // *** this is doing nothing now
    
        if (al.GetTagType(tName, tagTypeName))  // reads w/o bc return false
        {    
            if (tagTypeName == 'Z')             // verify that type is correct
            {               
                al.GetTag(tName, bc);           // extract the tag value from entry

                if (bc.compare(bc_prev) != 0)   // bc is a new barcode
                {     
                    if (nobc_flag)              // store # of anchored reads w/o barcodes
                    { 
                        loAnchored = anchoredR; 
                        nobc_flag = false;
                    } 
                    else                        // output # of anchored reads w/ barcodes
                    { 
                        file << anchoredR << "\n"; 
                    } 
                    anchoredR = 0;              // reset counts for this well
                }                       

                bc_prev = bc;  
            }   
            else 
            {
                printf("Warning: tagType is not Z - entry ignored \n");
            }  
        } 

        if (read.anchored)   // determine if the read is anchored 
        {   
            totalAnchored ++;
            anchoredR ++;
        }

        if (read.clipped)    // record clipped reads
        {
            clipcount ++;
        }

        if (nobc_flag)      // record through reads w/o bc
        {  
            loReads ++;
        }

        totalR ++;
    }

    file << anchoredR << "\n";  // save # of anchored reads for the last bc

    // map_to_file(bamFname, MapWell, loReads, loAnchored); 

    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout<<"Duration: "<< duration <<" s.\n";

    reader.Close();
    file.close();

    printf("* %u ", totalAnchored);
    printf("out of %u reads are anchored.\n", totalR);
    printf("* %u ", clipcount);
    printf("out of %u reads are clipped.\n", totalR);

    // To test how long: use
    //std::clock_t start;
    //double duration;
    //start = std::clock();
    // /*code here*/
    //duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    //std::cout<<"Duration: "<< duration <<" s.\n";
}

